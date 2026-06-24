library(atsim)

library(tidyverse)
library(sf)
library(rnaturalearth)


##keep the same seed for the survival and battery functions
set.seed(42)

##Create the fish survival function
## shape is less than one so more juvenile mortality
fish_wb_p = get_weibull_parm(14*365,1/300,upper_bound = c(0.9,5),parallel=FALSE)
fish_s = survival_function(pweibull,list(scale=fish_wb_p$scale,shape=fish_wb_p$shape))

batt_days = c(300,600,1500,2500)

batt_wb_p = lapply(batt_days,function(x){get_weibull_parm(x,0.95,lower_bound = c(1,0),upper_bound = c(10,10),parallel=FALSE)})
batt_s = lapply(batt_wb_p,function(x){survival_function(pweibull,list(scale=x$scale,shape=x$shape))})
names(batt_s) = batt_days

saveRDS(fish_wb_p,"fish_wb_p.rds")
saveRDS(batt_wb_p,"batt_wb_p.rds")

saveRDS(batt_s,"batt_s.rds")
saveRDS(fish_s,"fish_s.rds")

map = ne_countries(scale=10,returnclass = "sf")
sf::sf_use_s2(FALSE)
atlantic_bbox <- st_bbox(c(xmin=-80,
                           ymin=42,
                           xmax=-15,
                           ymax=62),crs=st_crs(map))
atlantic_bbox2 = st_transform(atlantic_bbox,crs="+proj=laea +lat_0=50.75 +lon_0=-48.88 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
bbox_polygon <- st_as_sfc(atlantic_bbox2)
map2 = st_transform(map,crs="+proj=laea +lat_0=50.75 +lon_0=-48.88 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

atlantic_map <- st_intersection(map2, bbox_polygon) |>
    st_combine()
sf::sf_use_s2(TRUE)
                         

map_filtered = st_transform(atlantic_map,crs="+proj=laea +lat_0=50.75 +lon_0=-48.88 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
saveRDS(map,"map.rds")
saveRDS(map_filtered,"map_filtered.rds")
