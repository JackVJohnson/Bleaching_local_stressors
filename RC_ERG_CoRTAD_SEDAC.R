#create wds
popwd = "C:/Users/Jack V Johnson/Desktop/Population"
Bleaching_data_directory="C:/Users/Jack V Johnson/Desktop/Bleaching_synergistic"

library(raster)
library(dismo)
library(rgdal)
library(dplyr)
library(tictoc)

#load raster

setwd(popwd)

Pop <- raster(paste0(getwd(), "/gpw_v4_population_count_rev11_2010_2pt5_min.tif", 
                     sep = ""))

#look at the info 
Pop

#define spatial geogrpahical projections

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#project for human pop data
proj4string(Pop) <- crs.geo

#get bleaching data

setwd(Bleaching_data_directory)

Bleaching_data <- read.csv(file="RC_ERG_SST_DHW.csv", header = TRUE, sep=",")
Bleaching_data <- select(Bleaching_data, -c(X.1,X))

#sites with coords
#site name lat lon
w.coord <- Bleaching_data[,c(3,5,6)]
head(w.coord)

#define coordiantes columns
coordinates(w.coord)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(w.coord) <- crs.geo
summary(w.coord)

tic()
# Extract human population density for the world within 100 km buffer
Bleaching_data$human100km <- extract(Pop, w.coord, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
Bleaching_data$human50km <- extract(Pop, w.coord, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
Bleaching_data$human25km <- extract(Pop, w.coord, buffer = 25000, fun = sum)
toc()

# merge the human population with bleaching + cortad data
write.csv(Bleaching_data, "RC_ERG_SST_DHW_Pop.csv")
