args <- commandArgs(trailingOnly = TRUE)
start_index = as.numeric(args[1])
seed = as.numeric(args[2])
ann_surv = as.numeric(args[3])

library(atsim)

library(tidyverse)
library(sf)
library(rnaturalearth)

##End time for the overall simulation
end_time = as.POSIXct("2025-12-31 23:59:59",tz="UTC")
release_datetime = as.POSIXct(sample(seq(as.POSIXct("2010-06-01 12:00:00","UTC"),as.POSIXct("2025-06-01 12:00:00",tz="UTC"),by="year"),1),tz="UTC")

##constant survival function
fish_s_gen <- function(survival_prob){
    fish_s <- function(t){
        survival_prob^t
    }
}
daily_surv = get_time_surv(365,ann_surv)
fish_s = fish_s_gen(daily_surv)

batt_wb_p = readRDS("batt_wb_p.rds")
batt_s = readRDS("batt_s.rds")


##Set the seed for everything else
set.seed(seed)

##Add battery info
batt_days = c(300,600,1500,2500)

batt_length = sample(batt_days,1,prob=c(0.10,0.20,0.30,0.5),replace=TRUE)
ReleaseDateTime = release_datetime

map = readRDS("map.rds")
map_filtered = readRDS("map_filtered.rds")


start_crs = st_crs(map)
mid_crs = st_crs(map_filtered)
end_crs = start_crs

rec = read_csv("./rec_lite.csv")
rec_sf1 = st_as_sf(rec,coords=c("ReceiverLongitude","ReceiverLatitude"),crs=start_crs)

## Logistic detection probability function in metres
pdrf <- function(dm, b = c(0.75, -1 / 300)) {
    p <- 1 / (1 + exp(-(b[1] + b[2] * dm)))
    return(p)
}

rec_info = list(rec_geometry=rec_sf1$geometry,rec_id=rec_sf1$Receiver,rec_deploy_time=rec_sf1$DeployedDateTime,rec_ret_time=rec_sf1$RetrievedDateTime,
                detection_fn=pdrf)

crs1 = st_crs(map)
crs2 = st_crs(map_filtered)
xlim = c(-56,-42)
ylim = c(43,61)
ylim2 = c(43,55)
offshore_b_gen <- function(polygon,dist,ylim,xlim,crs1,crs2,n,change_time){
    polygon2 = st_transform(polygon,crs=crs2)
    buffer_p = st_buffer(polygon2,dist)
    buffer_p1 = st_transform(buffer_p,crs=crs1)
    bbox_coords <- c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
    rectangle_sfc <- st_as_sfc(st_bbox(bbox_coords))
    st_crs(rectangle_sfc) = crs1
    rectangle_sfc2 = st_transform(rectangle_sfc,crs=crs2)
    diff_p = st_difference(rectangle_sfc2,buffer_p)
    
    samp = st_sample(diff_p,n+1)
    samp2 = st_transform(samp,crs=crs1)
    samp_c = st_coordinates(samp2)
    df1 = data.frame(x=samp_c[1:n,1],y=samp_c[1:n,2])
    df = as.data.frame(cbind(df1,change_time))
    names(df) = c("x","y","change_datetime")
    df

}

inshore_b_gen <- function(polygon,dist,ylim,xlim,crs1,crs2,n,change_time){
    polygon2 = st_transform(polygon,crs=crs2)
    buffer_p = st_buffer(polygon2,dist)
    diff_1 = st_difference(buffer_p,polygon)
    bbox_coords <- c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
    rectangle_sfc <- st_as_sfc(st_bbox(bbox_coords))
    st_crs(rectangle_sfc) = crs1
    rectangle_sfc2 = st_transform(rectangle_sfc,crs=crs2)
    diff_2 = st_difference(rectangle_sfc2,diff_1)


    diff_2 <- st_difference(rectangle_sfc2, diff_2)

    samp = st_sample(diff_2,n+1)
    samp2 = st_transform(samp,crs=crs1)
    samp_c = st_coordinates(samp2)
    df1 = data.frame(x=samp_c[1:n,1],y=samp_c[1:n,2])
    df = as.data.frame(cbind(df1,change_time))
    names(df) = c("x","y","change_datetime")
    df



}

migration_b_gen <- function(polygon,dist,ylim,xlim,crs1,crs2,n,change_time,day_offset=365.25/2){
    offshore_time = change_time+60*60*24*(day_offset)
    inshore = inshore_b_gen(polygon,dist,ylim,xlim,crs1,crs2,n,change_time)
    offshore = offshore_b_gen(polygon,dist,ylim,xlim,crs1,crs2,n,offshore_time)
    df = rbind(inshore,offshore)
    df
}

b_dates = seq(as.POSIXct("2010-01-01 12:00:00","UTC"),as.POSIXct("2026-01-01 12:00:00",tz="UTC"),by="year")

convert_latlon_to_m <- function(lat,lon,crs){
    point = st_as_sf(data.frame(y=lon,x=lat),coords=c("x","y"),crs=st_crs(ne_states()))
    point2 = st_transform(point,crs=crs)
    st_coordinates(point2)[1,]
}

library(sf)



age = NA
style = NA


batt_fn = batt_s[[as.character(batt_length)]]
rel_date = ReleaseDateTime
##Age is in days
age = runif(1,4*365,12*365)
style = sample(c("inshore","offshore","migration"),1)
if(style == "inshore"){
    cur_b = inshore_b_gen(map_filtered,30000,ylim2,xlim,crs1,crs2,length(b_dates),b_dates)
}else if(style == "offshore"){
    cur_b = offshore_b_gen(map_filtered,30000,ylim,xlim,crs1,crs2,length(b_dates),b_dates)
}else{
    cur_b =migration_b_gen(map_filtered,30000,ylim2,xlim,crs1,crs2,length(b_dates),b_dates)
}
a_mat = matrix(c(0.00001,0.00002,0.000001,0.00002),nrow=2,ncol=2)
##a_mat = matrix(c(runif(1,0.0000001,0.00001),runif(1,-0.00001,0.00001),runif(1,-0.00001,0.00001),c(runif(1,0.0000001,0.00001))),nrow=2,ncol=2)
cur_opts = list(a=a_mat)
curr_indi = NULL

##some release points clip into land
##init_coord = sf::st_as_sf(data.frame(x = rel_loc[1], y = rel_loc[2]), 
                                        # #                         coords = c("x", "y"), crs = st_crs(map))
##if (any(sf::st_within(init_coord, map, sparse = FALSE)) == 
##    TRUE) {
## Give a random release location
rel_u = runif(1)
if(rel_u < 0.5){
    rel_locY = inshore_b_gen(map_filtered,30000,ylim2,xlim,crs1,crs2,1,b_dates[1])
}else{
    rel_locY = offshore_b_gen(map_filtered,30000,ylim2,xlim,crs1,crs2,1,b_dates[1])
}
rel_loc = c(rel_locY$x,rel_locY$y)
##}

                                        # Loop continues as long as 'my_object' is NULL or carries an error class
while (is.null(curr_indi) || inherits(curr_indi, "try-error")) {
    
                                        # try() captures errors without stopping the script
    curr_indi<- try({
        simulate_individual(start_index,rel_date,rel_loc,age,fish_s,batt_fn,rec_info,final_rast,map_filtered,cur_b,start_crs=start_crs,mid_crs=mid_crs,end_crs=end_crs,ou_opts=cur_opts,timestep_per_day=0.25)
    }, silent = FALSE) # Suppresses error messages in the console
}
if(!dir.exists("./simstoreC/")){
    dir.create("./simstoreC/")
}
curr_indi$fish_info$batt_length = batt_length
curr_indi$fish_info$age = age
saveRDS(curr_indi,paste0("./simstoreC2/indi_",ann_surv,"_",str_pad(start_index,5,pad="0"),".rds"))



