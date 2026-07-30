
#' Calculate the survival needed to hit an annual survival target
#'
#' @param num_time_periods the number of time periods (e.g, 365 for days)
#' @param annual_survival the annual survival
#'
#' @export
get_time_surv <- function(num_time_periods,annual_survival){
    exp(log(annual_survival)/num_time_periods)
}

#' Generate the weibull hazard function for days with a given shape/scale
#' @param shape the shape parameter
#' @param scale the scale parameter
weibull_h_gen <- function(shape,scale){

    h <- function(t){
        dweibull(t/1000,shape,scale)/(1-pweibull(t/1000,shape,scale))
    }
    h
}

#' Find Weibull distribution parameters for battery survival
#'
#'  Assumes that shape must > 1
#' 
#' @param days number of days for that survival percent to be matched
#' @param surv_per the survival percentage to hit
#' @param init_parm the initial parameters to use
#' 
#' @export
get_weibull_parm <- function(days,surv_per,init_parm=c(1,1),lower_bound=c(0,0),upper_bound=c(5,5),maxiter=1000,numIslands=8,parallel=TRUE){
    cdf_per = 1-surv_per

    wp_opt <- function(parm){
        shape = parm[1]
        scale = parm[2]
        
        pwe = pweibull(days/1000,shape,scale)

        -(cdf_per-pwe)^2
    }

        GAopt = GA::gaisl(type="real-valued",fitness=wp_opt,lower=lower_bound,
                      upper=upper_bound,maxiter=maxiter,popSize=500,
                      numIslands = numIslands,migrationRate = 0.2,parallel = parallel)
    sol = GAopt@solution[1,]
    ret = list(shape=sol[1],scale=sol[2],opt=GAopt)
    ret
}

#' Vitality (2009 Li & Anderson) F dist
#'
#' One minus the survival equation described
#'
#' @param q the quantile (time)
#' @param r wear rate
#' @param s standard deviation in wear rate
#' @param k rate of accidental failure
#' @param u standard deviation in accidental failure
#' 
#' @export
pvitality <- function(q,r,s,k,u){
    raw_surv <- function(t){
        ncdf1 = (1-r*t)/sqrt(u^2+s^2*t)
        ncdf2 = -((2*u^2*r)/(s^2)+r*t+1)/sqrt(u^2+s^2*t)
        e1 = (2*u^2*r^2)/s^4 + (2*r)/s^2
        ret = (pnorm(ncdf1)-exp(e1)*pnorm(ncdf2))*exp(-k*t)
        ret
    }

    ##Normalizing so it always starts at 1
    time_zero_surv = raw_surv(0)

    ret = raw_surv(q)/time_zero_surv
    
    
    1-ret 
}

#' Quantile of Vitality 2009
#'
#' @param p the target percentile
#' @param r wear rate
#' @param s standard deviation in wear rate
#' @param k rate of accidental failure
#' @param u standard deviation in accidental failure
#'
#' @export
qvitality <- function(p,r,s,k,u){

    qp <- function(q,pp){
        pp-pvitality(q,r,s,k,u)
    }

    qq = sapply(p,function(x){
        uniroot(qp,c(0,100000),pp=x)$root
    })

    qq
}


#' Random values of Vitality 2009
#'
#' @param n the number required
#' @param r wear rate
#' @param s standard deviation in wear rate
#' @param k rate of accidental failure
#' @param u standard deviation in accidental failure
#'
#' @export
rvitality <- function(n,r,s,k,u){
    unival = runif(n)
    outs = qvitality(unival,r,s,k,u)
    outs
}

#' Find Vitality 2009 parameters for battery survival
#'
#' @param days1 number of days for the first survival percent to be matched
#' @param surv_per1 the first survival percent to be matched
#' @param maxiter the number of iterations for the genetic algoritim
#' @param numIslands the number of islands for the genetic algorithm to find the values
#' @param upper_bound the upper bound for the values of r, s, k and u
#' 
#' @export
get_vitality_parm <- function(days1,surv_per1,maxiter=1000,numIslands=8,upper_bound=c(5,5,5,0.05)){

    cdf_per1 = 1-surv_per1

    
    opt_fn <- function(parm){
        e1 = (cdf_per1-pvitality(days1/1000,parm[1],parm[2],parm[3],parm[4]))
        error = -((e1)^2)
        error
    }

    GAopt = GA::gaisl(type="real-valued",fitness=opt_fn,lower=c(0,0,0,0),
                      upper=upper_bound,maxiter=maxiter,popSize=500,
                      numIslands = numIslands,migrationRate = 0.2)
    sol = GAopt@solution[1,]
    ret = list(r=sol[1],s=sol[2],k=sol[3],u=sol[4],opt=GAopt)

    ret
}
               

#' Get the survival given some set of parameters
#' @param F_dist the F distribution to get survival function of
#' @param parameters parameters of the distribution as a list
#' @param t_adj the amount to adjust t 
#'
#' @export
survival_function <- function(F_dist,parameters,t_adj=1000){
    s_f <- function(t){
        args = c(list(q=t/t_adj),parameters)
        1-do.call(F_dist,args)
    }
}



#' Simulate a 2D Ornstein-Uhlenbeck process with boundary
#'
#' This merges elements of glatos' crw_in_polygon and adeHabitatLT simm.mou. You can simulate
#' points corresponding to a Ornstein-Uhlenbeck process with an attraction point but that won't cross
#' polygon boundries. The boundries are avoided simply by trying to find another point that doesn't cause the path
#' to cross the boundry, so in some instances it can have trouble and will fail if it can't find a good point.
#' The attraction points can be set seperate for every time step potentially and so can move.
#' 
#' @param polygon the sf polygon boundry the movement should not cross (polygon crs should be in metres)
#' @param timesteps the timesteps of the movement
#' @param b vector or matrix of attraction point(s), if matrix each row is the attraction point at the timestep 
#' @param a 2x2 matrix controlling the force of the attraction
#' @param x0 the intial start point
#' @param sigma the 2x2 matrix controlling the random noise in the movement
#' @param speed max speed in km/h
#' @param maxtries the number of tries to keep finding a point outside the boundary
#'
#' @export
sim_ou_boundry <- function(polygon,timesteps,b,a=diag(0.5,2),x0,sigma=matrix(c(50,0,50,0),nrow=2),speed=2,maxtries=100){

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
    buffered_land = sf::st_buffer(polygon,50)
    buffered_land = sf::st_union(buffered_land)
    
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
        zl_coord = st_as_sf(data.frame(x=zlast[1],y=zlast[2]),coords=c("x","y"),crs=st_crs(polygon))
        c_path = rbind(zl_coord,zcoord) |>
            st_combine() |>
            st_cast("LINESTRING")
        check_ppaths <- function(point_a,point_b){
            c_path = rbind(point_a,point_b) |>
                st_combine() |>
                st_cast("LINESTRING")
            any(sf::st_crosses(c_path,polygon,sparse=FALSE))
        }

        t_dist = max_dist
        if(any(sf::st_crosses(c_path,polygon,sparse=FALSE)) == TRUE){
            ##instead get a random point from the
            b_coord = st_as_sf(data.frame(x=b[1],y=b[2]),coords=c("x","y"),crs=st_crs(polygon))
            ntry = 0
            
            repeat{
                zl_buff = sf::st_buffer(zl_coord,t_dist)
                zl_diff = sf::st_difference(zl_buff,buffered_land)
        
                zsamp = sf::st_sample(zl_diff,100)
               
            
            
                dists = as.numeric(sf::st_distance(b_coord,zsamp)[1,])
                weights = 1/(dists)
                zsamp_coords = sf::st_coordinates(zsamp)
                zsampg = zsamp_coords[sample(1:nrow(zsamp_coords),1,prob=weights),]
                zcoord = st_as_sf(data.frame(x=zsampg[1],y=zsampg[2]),coords=c("x","y"),crs=st_crs(polygon))
                if(check_ppaths(zl_coord,zcoord) == FALSE){
                    break
                }
                ntry = ntry + 1
                t_dist = t_dist/2
                ##if things are reallllly bad
                if(ntry > 5){
                    zcoord = zl_coord
                    break
                }
            }

            }
        
        
 
        ret = list(zcoord = zcoord)
        #zcoord
        ret
    }
    
    Z = init_coord
    pb <- txtProgressBar(min = 1, n, style = 3)
    retret = list()
    for(i in 2:n){
        setTxtProgressBar(pb, i)
        Zt = Z[i-1,]
        dt = tdiff[i-1]
        dz = data.frame(x=rnorm(1,sd=sqrt(dt)),y=rnorm(1,sd=sqrt(dt)))
        W = as.vector(as.matrix(dz)%*%sigma)
        Zttt = make_z(a,as.vector(b[i,]),as.vector(st_coordinates(Zt)),dt,W,ms_speed)
        Zt = Zttt$zcoord
        retret[[i-1]] = Zttt
        Z =rbind(Z,Zt)
    }
    close(pb)


    Z$timesteps = timesteps
    Z
    ##list(Z=Z)
}

#' Add transmissions to a path that has known timesteps at a location
#'
#' @param path the sf points to add transmissions
#' @param timesteps the POSIXct timesteps of the path
#' @param delayRNG the min/max delay time from transmission to the next
#' @param burstDur number of seconds for the duration of the transmission
#' @param tz timezone of the time steps
path_transmissions <- function(path,timesteps,delayRNG = c(60,180),burstDur=5,tz="UTC"){

    sec_diff = as.numeric(difftime(timesteps[length(timesteps)],timesteps[1],units="secs",tz=tz))
    max_tran = sec_diff/(delayRNG[1]+burstDur)
    trans_t = runif(max_tran,delayRNG[1],delayRNG[2])+burstDur
    ctrans_t = cumsum(trans_t)
    ntimesteps = as.POSIXct(timesteps[1] + ctrans_t,tz=tz)
    ntimesteps = ntimesteps[ntimesteps <= timesteps[length(timesteps)]]
    xxx = approx(as.numeric(timesteps),sf::st_coordinates(path)[,"X"],xout=as.numeric(ntimesteps))$y
    yyy = approx(as.numeric(timesteps),sf::st_coordinates(path)[,"Y"],xout=as.numeric(ntimesteps))$y
    ret = data.frame(x=xxx,y=yyy,transmissionTime=ntimesteps)
    ret_sf = sf::st_as_sf(ret,coords=c("x","y"),crs=sf::st_crs(path))
    ret_sf 
}



#' Determine which transmissions get detected on a path
#'
#' @param path the sf geometry from path_transmissions containing the points at each transmission time
#' @param transmissionTime the vector of transmission times from path_transmissions output
#' @param rec_geometry the sf geometry of points for receivers
#' @param rec_id the vector of receiver ids (e.g., Receiver)
#' @param rec_deploy_time the vector of receiver deployment times (e.g., DeployedDateTime)
#' @param rec_ret_time the vector of receiver retrieval times (e.g., RetrivedDateTime)
#' @param detection_fn a function giving the probability of a transmission being detected in metres
#' @param fish_id the id you want to add to the detections
#' @param end_crs (optional) the crs to use at the end so that you can transform the receiver and fish locations on output
#' 
#' @export
#' 
detections_on_path <- function(path,transmissionTime,rec_geometry,rec_id,rec_deploy_time,rec_ret_time,detection_fn,fish_id,end_crs=NULL){

    ##Find the detection range distance that has probability be essentially zero
    zero_dist = uniroot(function(x){detection_fn(x)-1e-18},c(0,50000))$root

    ##Buffer the receivers and then find points in that buffer
    ##path_buffer = sf::st_buffer(path,zero_dist)
    rec_buffer = sf::st_buffer(rec_geometry,zero_dist)
    overlap_b = sf::st_intersects(path,rec_buffer)
    overlappers = as.data.frame(overlap_b)
    dists = sf::st_distance(path[overlappers[,1],],rec_geometry[overlappers[,2],],by_element = TRUE)
    units(dists) = units::make_units(m)
    pdists = detection_fn(as.numeric(dists))
    suc_trans = rbinom(length(pdists),1,pdists)
    pos_trans = overlappers[which(suc_trans == 1),]
    if(nrow(pos_trans) > 0){

        rec_geometryL = rec_geometry[pos_trans[,2],]
        pathL = path[pos_trans[,1],]
        if(!is.null(end_crs)){
            rec_geometryL = sf::st_transform(rec_geometryL,crs=end_crs)
            pathL = sf::st_transform(pathL,crs=end_crs)
        }

    ##Check that the transmission time is within the receivers deployment range
        df_time = data.frame(DetectionDateTime=transmissionTime[pos_trans[,1]],DeployedDateTime=rec_deploy_time[pos_trans[,2]],RetrievedDateTime=rec_ret_time[pos_trans[,2]],Receiver=rec_id[pos_trans[,2]],rec_x=sf::st_coordinates(rec_geometryL)[,"X"],rec_y=sf::st_coordinates(rec_geometryL)[,"Y"],fish_x=sf::st_coordinates(pathL)[,"X"],fish_y=sf::st_coordinates(pathL)[,"Y"],fish_id=fish_id)
        df_time = dplyr::filter(df_time,DetectionDateTime >= DeployedDateTime & DetectionDateTime <= RetrievedDateTime)

    }else{
                df_time = data.frame(DetectionDateTime=POSIXct(0),DeployedDateTime=POSIXct(0),RetrievedDateTime=POSIXct(0),Receiver=character(0),rec_x=numeric(0),rec_y=numeric(0),fish_x=numeric(0),fish_y=numeric(0),fish_id=character(0))

    }
    df_time
}

#' Simulate one individual's detections and track.
#'
#' @param fish_id the ID to give the individual
#' @param release_datetime the release date time of the individual (POSIXct)
#' @param release_location vector of the release location of the individual
#' @param start_age the age at release (in days)
#' @param individual_s_fn the survival function of the individual (in days)
#' @param battery_s_fn the survival function of the transmitter battery (in days)
#' @param rec_info named list containing the rec_geometry, rec_id, rec_deploy_time, rec_ret_time, and detection_fn for detections_on_path
#' @param polygon the polygon for sim_ou_boundary (must have mid_crs)
#' @param b_frame the attrator points in the start_crs, data frame with x,y, and change_datetime for when that attrator point applies 
#' @param ou_opts optional named list giving the parameters of sim_ou_boundry other than b and timesteps
#' @param start_crs the crs of the release_location 
#' @param mid_crs the crs used for sim_ou_boundry
#' @param end_crs the final output crs
#' @param end_time the end date time of the simulation
#' @param delayRNG the min/max delay time from transmission to the next
#' @param burstDur number of seconds for the duration of the transmission
#' @param timestep_per_day the number of tries for sim_ou_boundry
#' @export
simulate_individual <- function(fish_id,release_datetime,release_location,start_age,individual_s_fn,battery_s_fn,rec_info,raster,polygon,b_frame,
                                ou_opts=NULL,start_crs,mid_crs,end_crs,end_time=as.POSIXct("2025-12-31 23:59:59",tz="UTC"),
                                delayRNG=c(60,180),burstDur=5,timestep_per_day=0.75){


    ##Calculate the death date of an individual
    fish_up_s <- function(q){individual_s_fn(start_age+q)/individual_s_fn(start_age)}
    fish_u = runif(1)
    fish_s_fn = function(x){fish_up_s(x)-fish_u}
    fish_s_opt = uniroot(function(x){fish_up_s(x)-fish_u},c(0,20*365))$root
    death_day = as.POSIXct(as.numeric(release_datetime)+fish_s_opt*24*60*60,tz="UTC")

    ##calculate the failure day of the transmitter battery
    batt_u = runif(1)
    batt_s_opt = uniroot(function(x){battery_s_fn(x)-batt_u},c(0,20*365))$root
    batt_failure_day =  as.POSIXct(as.numeric(release_datetime)+batt_s_opt*24*60*60,tz="UTC")

    ##Create the b matrix for sim_ou_boundry
    max_day = min(death_day,batt_failure_day,end_time)
    start_to_max = difftime(max_day,release_datetime,units="days")
    timesteps = sort(as.POSIXct(sample(release_datetime:max_day,floor(as.numeric(start_to_max)*timestep_per_day)),tz="UTC"))
    b_f_sf = sf::st_as_sf(b_frame,coords=c("x","y"),crs=start_crs)
    b_f_sf = sf::st_transform(b_f_sf,crs=mid_crs)
    t_df = data.frame(timestep=timesteps)
    b_join = dplyr::join_by(closest(timestep >= change_datetime))
    b_df = dplyr::left_join(t_df,b_f_sf,by=b_join)
    b_mat = sf::st_coordinates(b_df$geometry)

    x0_df = data.frame(x=release_location[1],y=release_location[2])
    x0_sf = sf::st_as_sf(x0_df,coords=c("x","y"),crs=start_crs)
    x0_sf = sf::st_transform(x0_sf,crs=mid_crs)
    x0 = sf::st_coordinates(x0_sf)

    if(is.null(ou_opts)){
        ou_opts = list()
    }

    if(is.null(ou_opts$a)){
        ou_opts$a = matrix(runif(4,max=0.0001),nrow=2,ncol=2)
    }

    if(is.null(ou_opts$sigma)){
        ou_opts$sigma = matrix(runif(4,max=500),nrow=2,ncol=2)
    }

    if(is.null(ou_opts$speed)){
        ou_opts$speed = 0.75
    }

    ou_process = sim_ou_boundry(polygon,timesteps,b_mat,ou_opts$a,x0[1,],ou_opts$sigma,ou_opts$speed)
    
    ##add more points if the battery failure day is greater than the death day
    if(batt_failure_day > death_day){
        new_row = sf::st_as_sf(data.frame(geometry=ou_process$geometry[nrow(ou_process)],timesteps=batt_failure_day),crs=sf::st_crs(ou_process))
        ou_process2 = rbind(ou_process,new_row)
    }else{
        ou_process2 = ou_process
    }

    transd_path = path_transmissions(ou_process2$geometry,ou_process2$timesteps,delayRNG,burstDur,tz="UTC")
    
    
    rec_df = data.frame(rec_id=rec_info$rec_id,rec_deploy_time=rec_info$rec_deploy_time,rec_ret_time=rec_info$rec_ret_time,rec_row=1:length(rec_info$rec_id))
   
    rec_geo = rec_info$rec_geometry
    rec_geo = sf::st_transform(rec_geo,crs=mid_crs)
    
    
    print("Finding detections...")
    detections = detections_on_path(transd_path$geometry,transd_path$transmissionTime,rec_geo,rec_df$rec_id,rec_df$rec_deploy_time,rec_df$rec_ret_time,rec_info$detection_fn,fish_id,end_crs=end_crs)

     x0_sf_end = sf::st_transform(x0_sf,crs=end_crs)

    x0_end = sf::st_coordinates(x0_sf_end)[1,]

    fish_info = data.frame(fish_id=fish_id,ReleaseDateTime=release_datetime,DeathDateTime=death_day,BatteryFailDate=batt_failure_day,rel_x=x0_end[1],rel_y=x0_end[2],delay_min=delayRNG[1],delay_max=delayRNG[2],burstDur=burstDur)

    out_path = sf::st_transform(ou_process,crs=end_crs)
    b_df_out = sf::st_transform(sf::st_as_sf(b_df),crs=end_crs)
    ret = list(fish_info=fish_info,detections=detections,path=out_path,b_df=b_df_out,ou_opts=ou_opts)
    ret
    
}

#' Plot an indivdual's path and attraction points
#'
#' @param individual output of simulate_individual
#' @export
plot_individual <- function(individual,ylim=NULL,xlim=NULL){

    map = rnaturalearth::ne_coastline(scale=10)

    path = sf::st_transform(individual$path,crs=sf::st_crs(map))
    b_df = individual$b_df

    path_coords = as.data.frame(sf::st_coordinates(path))
    b_coords = as.data.frame(sf::st_coordinates(sf::st_transform(b_df$geometry,crs=sf::st_crs(map))))

    if(is.null(ylim)){
        ylim =c(min(path_coords$Y),max(path_coords$Y))
    }

    if(is.null(xlim)){
        xlim =c(min(path_coords$X),max(path_coords$X))
    }

    
    pl = ggplot2::ggplot() + geom_sf(data=map)  + geom_path(aes(color=path$timesteps,x=X,y=Y),data=path_coords) +  geom_point(data=path_coords[1,],aes(x=X,y=Y),shape=23,color="orange",fill="orange",size=2) + geom_point(data=path_coords[nrow(path_coords),],aes(x=X,y=Y),shape=23,color="green",fill="green",size=2) + geom_point(data=b_coords,aes(x=X,y=Y),color="red",shape=24)
    if(nrow(individual$detections) > 0){
        qd = individual$detections |>
            dplyr::group_by(rec_x,rec_y) |>
            dplyr::summarise(n_det=n())
        
       pl = pl + geom_point(data=qd,aes(x=rec_x,y=rec_y,size=sqrt(n_det)))
    }
    pl = pl + ggplot2::coord_sf(xlim=xlim,ylim=ylim,crs=sf::st_crs(map))

    pl

}

#' Calculate the emperical survival function from the simulation
#'
#' @param ReleaseDateTime the datetime of each released individual
#' @param DeathDateTime the datetime of each individual's death
#'
#' @export
sim_emp_surv_fn <- function(ReleaseDateTime,DeathDateTime){

    n_ind = length(ReleaseDateTime)
    days_surv = difftime(DeathDateTime,ReleaseDateTime)
    ret_fn <- function(t,unit="days"){
        t = as.difftime(t,units=unit)
        
        sapply(t,function(x){sum(days_surv >= x)/n_ind})
    }


}

#' Build dataset from individual simulations
#'
#' Merge together individuals into one dataset after being
#' randomly sampled. Returns the combined fish_info dataframe
#' and combined detections of the sampled individuals.
#'
#' @param indiv_list list of objects from simulate_individuals
#' @param n the number of individuals to include in the data set
#' @param crs if not null the CRS to convert x and y coordinates to (assumes lat and lon from start)
#' 
#' @export
#'
build_sim_dataset <- function(indiv_list,n,crs=NULL){
    to_get = sample(1:length(indiv_list),n)

    got = indiv_list[to_get]
    fish_info = lapply(got,function(x){
        y = x$fish_info
        y$n_det = nrow(x$detections)
        y$b_fail_before_d = y$BatteryFailDate <= y$DeathDateTime
        y
    })
    fish_info <- data.table::rbindlist(fish_info, use.names = TRUE, fill = TRUE)
    if(is.numeric(fish_info$fish_id)){
        fish_info$fish_id = as.character(fish_info$fish_id)
    }
    
    detections = lapply(got,function(x){
        x$detections
    })
    detections <- data.table::rbindlist(detections, use.names = TRUE, fill = TRUE)

    if(!is.null(crs)){
        dsf1 = sf::st_as_sf(detections,coords=c("rec_x","rec_y"),crs=4326)
        dsf1t = sf::st_transform(dsf1,crs)
        d1c = sf::st_coordinates(dsf1t)
        detections$rec_x = d1c[,1]
        detections$rec_y = d1c[,2]

        dsf2 = sf::st_as_sf(detections,coords=c("fish_x","fish_y"),crs=4326)
        dsf2t = sf::st_transform(dsf2,crs)
        d2c = sf::st_coordinates(dsf2t)
        detections$fish_x = d2c[,1]
        detections$fish_y = d2c[,2]

        fsf = sf::st_as_sf(fish_info,coords=c("rel_x","rel_y"),crs=4326)
        fsft = sf::st_transform(fsf,crs)
        fc = sf::st_coordinates(fsft)
        fish_info$rel_x = fc[,1]
        fish_info$rel_y = fc[,2]
    }
        

    ret = list(fish_info=fish_info,detections=detections)
    ret
    
    
}
