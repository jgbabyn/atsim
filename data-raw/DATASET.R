## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
library(rnaturalearth)
library(sf)
library(tidyverse)

 na_map <- ne_states(country = c("united states of america", "canada","greenland","saint pierre and miquelon"), returnclass = "sf") 
  # 2. Filter out Hawaii and Alaska (if desired)
  # Use 'name' or 'name_en' to identify the regions
  map_filtered <- na_map |>
   filter(!(name %in% c("Hawaii","Alaska", "Puerto Rico", "Guam", 
                        "American Samoa", "United States Virgin Islands", 
                        "Northern Mariana Islands"))) |>
      st_make_valid() |>
      st_union()

map_filtered = st_transform(map_filtered,crs="+proj=laea +lat_0=50.75 +lon_0=-48.88 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
points = data.frame(x=c(-56,-52,-60),y=c(48,48,52))

points_sf = st_as_sf(points,coords=c("x","y"),crs=st_crs(ne_states()))
points_sf2 = st_transform(points_sf,crs=st_crs(map_filtered))
any(st_within(points_sf2,map_filtered,sparse=FALSE))
timesteps = seq(as.POSIXct("2010-01-01 12:00:00",tz="UTC"),as.POSIXct("2012-02-01 12:00:00",tz="UTC"),by="day")

convert_latlon_to_m <- function(lat,lon,crs){
    point = st_as_sf(data.frame(y=lon,x=lat),coords=c("x","y"),crs=st_crs(ne_states()))
    point2 = st_transform(point,crs=crs)
    st_coordinates(point2)[1,]
}

b = convert_latlon_to_m(-51,49,st_crs(map_filtered))
x0=convert_latlon_to_m(-52,47,st_crs(map_filtered))
in_point = point_to_sf(-54,48,st_crs(map_filtered))
out_point = point_to_sf(-52,47,st_crs(map_filtered))

point_to_sf <- function(lat,lon,crs){
    point = st_as_sf(data.frame(y=lon,x=lat),coords=c("x","y"),crs=st_crs(ne_states()))
    point2 = st_transform(point,crs=crs)
    point2
}

blong = matrix(b,nrow=length(timesteps),ncol=2,byrow=TRUE)
for(i in 1000:2000){
    blong[i,] = convert_latlon_to_m(-52,47,st_crs(map_filtered))
}


sim_ou_boundry2 <- function(polygon,timesteps,b,a=diag(0.5,2),x0,sigma=matrix(c(50,0,50,0),nrow=2),speed=2,maxtries=100){

    pg_unit = sf::st_crs(polygon, parameters = TRUE)$units_gdal
    if(pg_unit == "metre"){
        ms_speed = (speed*1000)/3600
        
    }else if(pg_unit != "metre"){
        stop("Please use a crs with metres")
    }

    if(is.vector(b)){
        b = matrix(b,nrow=length(timesteps),ncol=2,byrow=TRUE)
    }
    if(nrow(b) != length(timesteps)){
        stop("b does not have the correct number of rows")
    }
    init_coord = sf::st_as_sf(data.frame(x=x0[1],y=x0[2]),coords=c("x","y"),crs=st_crs(polygon))
    b_coords = sf::st_as_sf(data.frame(x=b[,1],y=b[,2]),coords=c("x","y"),crs=st_crs(polygon))
    if(any(sf::st_within(init_coord,polygon,sparse=FALSE)) == TRUE){
        stop("x0 is within the boundary")
    }

    if(any(sf::st_within(b_coords,polygon,sparse=FALSE)) == TRUE){
        stop("some value of b is within the boundary")
    }
    
    if(!all.equal(dim(a),c(2,2))){
        stop("a must be 2x2 matrix")
    }
    if(!all.equal(dim(sigma),c(2,2))){
        stop("sigma must be 2x2 matrix")
    }

   
    

    

    n = length(timesteps)
    tdiff = as.numeric(diff(timesteps,units="secs"))


    make_z <- function(a,b,zlast,dt,wlast,max_speed=ms_speed){
        z_raw1 = a%*%(b-zlast)*dt+wlast
        z_noise = wlast/dt
        z_dist = sqrt(z_raw1[1,1]^2+z_raw1[2,1]^2)
        max_dist = max_speed*dt
        if(z_dist > max_dist){
            z_raw = z_raw1*(max_dist/z_dist)
        }else{
            z_raw = z_raw1
        }
        point = zlast + z_raw
        zcoord = st_as_sf(data.frame(x=point[1,1],y=point[2,1]),coords=c("x","y"),crs=st_crs(polygon))
        dp = as.numeric(sf::st_distance(zcoord,polygon)[1])
        if(dp == 0){
            point = zlast - z_raw
            zcoord = st_as_sf(data.frame(x=point[1,1],y=point[2,1]),coords=c("x","y"),crs=st_crs(polygon))
        }
        ret = list(zcoord=zcoord,point=point,z_raw=z_raw,zback=z_raw*dt+wlast,wlast=wlast,dt=dt,a=a,b=b,zlast=zlast,z_dist=z_dist,max_dist=max_dist,z_raw1=z_raw1,dp=dp)
        #zcoord
        ret
    }
    
    Z = init_coord
    pb <- txtProgressBar(min = 1, n, style = 3)
    retret = list()
    dps = list()
    for(i in 2:n){
        setTxtProgressBar(pb, i)
        Zt = Z[i-1,]
        dt = tdiff[i-1]
        ntries = 0
        while(TRUE){
            dz = data.frame(x=rnorm(1,sd=sqrt(dt)),y=rnorm(1,sd=sqrt(dt)))
            W = as.vector(as.matrix(dz)%*%sigma)
            Zttt = make_z(a,as.vector(b[i,]),as.vector(st_coordinates(Zt)),dt,W,ms_speed)
            Zt = Zttt$zcoord
            dps[[i]] = Zttt$dp
            c_path = rbind(Z,Zt) |>
                st_combine() |>
                st_cast("LINESTRING")
            if(any(sf::st_crosses(c_path,polygon,sparse=FALSE)) == FALSE){
                retret[[i-1]] = Zttt
                break
            }
            ntries = ntries + 1
            #cat("On",as.character(i),"on the",as.character(ntries),"try.")
            if(ntries > maxtries){
                stop("stuck somewhere?")
            }
        }
        Z =rbind(Z,Zt)
    }
    close(pb)


    Z$timesteps = timesteps
    Z
    list(Z=Z,dps)
}


system.time(goof <- sim_ou_boundry2(map_filtered,a=matrix(c(0.00001,0.00002,0.000001,0.00002),nrow=2),timesteps,b=blong,x0=convert_latlon_to_m(-52,47,st_crs(map_filtered)),speed=0.75,sigma=matrix(c(200,0,200,0),nrow=2),maxtries=1))
goofls = goof |>
    arrange(timesteps) |>
    st_combine() |>
    st_cast("LINESTRING")

library(ggplot2)
init_point = point_to_sf(-52,47,st_crs(map_filtered))
attract_point = point_to_sf(-51,49,st_crs(map_filtered))
attract_point2 = point_to_sf(-54,48,st_crs(map_filtered))

zl_point = init_point
end_sf = attract_point2

minp = convert_latlon_to_m(-60,45,st_crs(map_filtered))
maxp = convert_latlon_to_m(-47,55,st_crs(map_filtered))
gpl = ggplot() + geom_sf(data=map_filtered) + geom_sf(data=goofls) + geom_sf(data=init_point,color="red") + geom_sf(data=attract_point,color="blue") + coord_sf(xlim=c(minp[1],maxp[1]),ylim=c(minp[2],maxp[2]),crs=st_crs(map_filtered)) 

rec = readRDS("/home/jbabyn/NCAT/TelemetryR_26/NCATexplorer/data/all_receivers.rds")

rec_sf = rec |>
    st_as_sf(coords=c("ReceiverLongitude","ReceiverLatitude"),crs=st_crs(ne_states())) |>
    st_transform(crs=st_crs(map_filtered)) |>
    select(geometry)

rec_cood = as.data.frame(st_coordinates(rec_sf))
names(rec_cood) = c("x","y")

     pdrf <- function(dm, b = c(0.5, -1 / 120)) {
       p <- 1 / (1 + exp(-(b[1] + b[2] * dm)))
       return(p)
     }
pdrf(c(50,100, 200, 300, 400, 500)) # view detection probs. at some distances

mytrans = path_transmissions(goof$Z,timesteps)
transmit = transmit_along_path(goof$Z,vel=0.75,sp_out=FALSE)
mydet = detect_transmissions(transmit,rec_cood,pdrf,sp_out=FALSE)

mypoly <- data.frame(
       x = c(0, 0, 1000, 1000),
       y = c(0, 1000, 1000, 0)
     )
     
     mypath <- crw_in_polygon(mypoly,
       stepLen = 100,
       nsteps = 50,
       sp_out = TRUE
     )
     
     plot(mypath, type = "l", xlim = c(0, 1000), ylim = c(0, 1000))
     
     # add receivers
     recs <- expand.grid(x = c(250, 750), y = c(250, 750))
     points(recs, pch = 15, col = "blue")
     
     # simulate tag transmissions
     mytrns <- transmit_along_path(mypath,
       vel = 2.0, delayRng = c(60, 180),
       burstDur = 5.0, sp_out = TRUE
     )
     points(mytrns, pch = 21) # add to plot
     
     # Define detection range function (to pass as detRngFun)
     # that returns detection probability for given distance
     # assume logistic form of detection range curve where
     #   dm = distance in meters
     #   b = intercept and slope
     pdrf <- function(dm, b = c(0.5, -1 / 120)) {
       p <- 1 / (1 + exp(-(b[1] + b[2] * dm)))
       return(p)
     }
     pdrf(c(100, 200, 300, 400, 500)) # view detection probs. at some distances
     
     # simulate detection
     mydtc <- detect_transmissions(
       trnsLoc = mytrns,
       recLoc = recs,
       detRngFun = pdrf,
       sp_out = FALSE
     )
     

fish_wb_p = get_weibull_parm(14*365,1/302,upper_bound=c(0.75,5))

fish_h = weibull_h_gen(fish_wb_p$shape,fish_wb_p$scale)

fish_s = survival_function(pweibull,list(scale=fish_wb_p$scale,shape=fish_wb_p$shape))

end_time = as.POSIXct("2025-12-31 23:59:59",tz="UTC")
release_datetime = as.POSIXct("2019-12-31 23:59:59",tz="UTC")

batt_wb_p = get_weibull_parm(3001,0.95,lower_bound = c(1,0),upper_bound = c(10,10))
battery_s_fn = survival_function(pweibull,list(scale=batt_wb_p$scale,shape=batt_wb_p$shape))

start_crs = st_crs(ne_states())
mid_crs = st_crs(map_filtered)
end_crs = start_crs
b_frame = data.frame(x=-51,y=49,change_datetime=as.POSIXct("2010-01-01 12:00:00"))
b_frame[2,] = c(-52,47,change_datetime = as.POSIXct("2020-04-01 12:00:00"))

release_location = c(-52,47)

polygon = map_filtered

burstDur = 5
delayRNG

rec_sf1 = st_as_sf(rec,coords=c("ReceiverLongitude","ReceiverLatitude"),crs=st_crs(ne_states()))

fish_id = 1
rec_info = list(rec_geometry=rec_sf1$geometry,rec_id=rec_sf1$Receiver,rec_deploy_time=rec_sf1$DeployedDateTime,rec_ret_time=rec_sf1$RetrievedDateTime,
                detection_fn=pdrf)


ou_opts = list(a = matrix(c(0.00001,0.00002,0.000001,0.00002),nrow=2,ncol=2),sigma = matrix(c(200,0,200,0),nrow=2,byrow=TRUE))
ou_opts = ind1$ou_opts
ou_opts$sigma = matrix(c(200,0,200,0),nrow=2,byrow=TRUE)
ind1 = simulate_individual(1,as.POSIXct("2019-04-05 12:53:53",tz="UTC"),c(-51,49),5*365,fish_s,battery_s_fn,rec_info,map_filtered,b_frame,
                           start_crs=start_crs,mid_crs=mid_crs,end_crs=end_crs,ou_opts=ou_opts)

ind2 = simulate_individual(1,as.POSIXct("2019-04-05 12:53:53",tz="UTC"),c(-51,49),5*365,fish_s,battery_s_fn,rec_info,map_filtered,b_frame,
                           start_crs=start_crs,mid_crs=mid_crs,end_crs=end_crs,ou_opts=ou_opts)

xlim = c(-60,-43)
ylim = c(46,58)

plot_individual(ind1,ylim,xlim)


ggplot() + geom_sf(data=map_filtered) + geom_sf(data=g)

zland = st_coordinates(st_line_sample(garp,1))[1,1:2]
ip = st_coordinates(init_point)[1,]
lp_norm <- function(x,p){
    sum(abs(x)^p)^(1/p)
}

grad_land = function(zlast,zland){
    suby = zland-zlast
    -suby/(lp_norm(suby,3))
}

fish_lite = fish |>
    mutate(batt_diff = difftime(est_bat_day,ReleaseDateTime,units="days")) |>
    select(ReleaseDateTime,ReleaseLatitude,ReleaseLongitude) |>
    filter(year(ReleaseDateTime) >= 2010)
write_csv(fish_lite,"./data-raw/fish_lite.csv")

rec = readRDS("/home/jbabyn/NCAT/TelemetryR_26/NCATexplorer/data/all_receivers.rds")

rec_l = rec[grep("^VR",rec$Receiver),]
write_csv(rec_l,"./data-raw/rec_lite.csv")
