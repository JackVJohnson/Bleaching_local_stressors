#create wds
# further population estimates 

popwd = "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/codes and files/worldpop"
Bleaching_data_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/codes and files"

library(raster)
library(dismo)
library(rgdal)
library(dplyr)
library(tictoc)

#load raster

setwd(popwd)


pop_2003 <- raster("ppp_2003_1km_Aggregated.tif")
pop_2004 <- raster("ppp_2004_1km_Aggregated.tif")
pop_2005 <- raster("ppp_2005_1km_Aggregated.tif")
pop_2006 <- raster("ppp_2006_1km_Aggregated.tif")
pop_2007 <- raster("ppp_2007_1km_Aggregated.tif")
pop_2008 <- raster("ppp_2008_1km_Aggregated.tif")
pop_2009 <- raster("ppp_2009_1km_Aggregated.tif")
pop_2010 <- raster("ppp_2010_1km_Aggregated.tif")
pop_2011 <- raster("ppp_2011_1km_Aggregated.tif")
pop_2012 <- raster("ppp_2012_1km_Aggregated.tif")
pop_2013 <- raster("ppp_2013_1km_Aggregated.tif")
pop_2014 <- raster("ppp_2014_1km_Aggregated.tif")
pop_2015 <- raster("ppp_2015_1km_Aggregated.tif")
pop_2016 <- raster("ppp_2016_1km_Aggregated.tif")
pop_2017 <- raster("ppp_2017_1km_Aggregated.tif")



#define spatial geogrpahical projections

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

#project for human pop data
proj4string(pop_2003) <- crs.geo

#get bleaching data

setwd(Bleaching_data_directory)

Bleaching_data <- read.csv(file="RC_ERG_SST_DHW_Pop_MPA.csv", header = TRUE, sep=",")
Bleaching_data <- select(Bleaching_data, -c(X.1,X))
summary(Bleaching_data$Year)

# subset data by year

df2003<- subset(Bleaching_data, Year == 2003)
df2004<- subset(Bleaching_data, Year == 2004)
df2005<- subset(Bleaching_data, Year == 2005)
df2006<- subset(Bleaching_data, Year == 2006)
df2007<- subset(Bleaching_data, Year == 2007)
df2008<- subset(Bleaching_data, Year == 2008)
df2009<- subset(Bleaching_data, Year == 2009)
df2010<- subset(Bleaching_data, Year == 2010)
df2011<- subset(Bleaching_data, Year == 2011)
df2012<- subset(Bleaching_data, Year == 2012)
df2013<- subset(Bleaching_data, Year == 2013)
df2014<- subset(Bleaching_data, Year == 2014)
df2015<- subset(Bleaching_data, Year == 2015)
df2016<- subset(Bleaching_data, Year == 2016)
df2017<- subset(Bleaching_data, Year == 2017)

# 2003

#sites with coords
#site name lat lon
coord2003 <- df2003[,c(3,5,6)]
head(coord2003)

#define coordiantes columns
coordinates(coord2003)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2003) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2003$pop100km <- extract(pop_2003, coord2003, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2003$pop50km <- extract(pop_2003, coord2003, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2003$pop25km <- extract(pop_2003, coord2003, buffer = 25000, fun = sum)
toc()

# 2004

#sites with coords
#site name lat lon
coord2004 <- df2004[,c(3,5,6)]
head(coord2004)

#define coordiantes columns
coordinates(coord2004)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2004) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2004$pop100km <- extract(pop_2004, coord2004, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2004$pop50km <- extract(pop_2004, coord2004, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2004$pop25km <- extract(pop_2004, coord2004, buffer = 25000, fun = sum)
toc()

# 2005

#sites with coords
#site name lat lon
coord2005 <- df2005[,c(3,5,6)]
head(coord2005)

#define coordiantes columns
coordinates(coord2005)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2005) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2005$pop100km <- extract(pop_2005, coord2005, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2005$pop50km <- extract(pop_2005, coord2005, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2005$pop25km <- extract(pop_2005, coord2005, buffer = 25000, fun = sum)
toc()

# 2006

#sites with coords
#site name lat lon
coord2006 <- df2006[,c(3,5,6)]
head(coord2006)

#define coordiantes columns
coordinates(coord2006)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2006) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2006$pop100km <- extract(pop_2006, coord2006, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2006$pop50km <- extract(pop_2006, coord2006, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2006$pop25km <- extract(pop_2006, coord2006, buffer = 25000, fun = sum)
toc()

# 2007

#sites with coords
#site name lat lon
coord2007 <- df2007[,c(3,5,6)]
head(coord2007)

#define coordiantes columns
coordinates(coord2007)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2007) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2007$pop100km <- extract(pop_2007, coord2007, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2007$pop50km <- extract(pop_2007, coord2007, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2007$pop25km <- extract(pop_2007, coord2007, buffer = 25000, fun = sum)
toc()

# 2008

#sites with coords
#site name lat lon
coord2008 <- df2008[,c(3,5,6)]
head(coord2008)

#define coordiantes columns
coordinates(coord2008)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2008) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2008$pop100km <- extract(pop_2008, coord2008, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2008$pop50km <- extract(pop_2008, coord2008, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2008$pop25km <- extract(pop_2008, coord2008, buffer = 25000, fun = sum)
toc()

# 2009

#sites with coords
#site name lat lon
coord2009 <- df2009[,c(3,5,6)]
head(coord2009)

#define coordiantes columns
coordinates(coord2009)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2009) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2009$pop100km <- extract(pop_2009, coord2009, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2009$pop50km <- extract(pop_2009, coord2009, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2009$pop25km <- extract(pop_2009, coord2009, buffer = 25000, fun = sum)
toc()

# 2010

#sites with coords
#site name lat lon
coord2010 <- df2010[,c(3,5,6)]
head(coord2010)

#define coordiantes columns
coordinates(coord2010)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2010) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2010$pop100km <- extract(pop_2010, coord2010, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2010$pop50km <- extract(pop_2010, coord2010, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2010$pop25km <- extract(pop_2010, coord2010, buffer = 25000, fun = sum)
toc()

# 2011

#sites with coords
#site name lat lon
coord2011 <- df2011[,c(3,5,6)]
head(coord2011)

#define coordiantes columns
coordinates(coord2011)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2011) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2011$pop100km <- extract(pop_2011, coord2011, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2011$pop50km <- extract(pop_2011, coord2011, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2011$pop25km <- extract(pop_2011, coord2011, buffer = 25000, fun = sum)
toc()


# 2012

#sites with coords
#site name lat lon
coord2012 <- df2012[,c(3,5,6)]
head(coord2012)

#define coordiantes columns
coordinates(coord2012)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2012) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2012$pop100km <- extract(pop_2012, coord2012, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2012$pop50km <- extract(pop_2012, coord2012, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2012$pop25km <- extract(pop_2012, coord2012, buffer = 25000, fun = sum)
toc()

# 2013

#sites with coords
#site name lat lon
coord2013 <- df2013[,c(3,5,6)]
head(coord2013)

#define coordiantes columns
coordinates(coord2013)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2013) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2013$pop100km <- extract(pop_2013, coord2013, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2013$pop50km <- extract(pop_2013, coord2013, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2013$pop25km <- extract(pop_2013, coord2013, buffer = 25000, fun = sum)
toc()

# 2014

#sites with coords
#site name lat lon
coord2014 <- df2014[,c(3,5,6)]
head(coord2014)

#define coordiantes columns
coordinates(coord2014)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2014) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2014$pop100km <- extract(pop_2014, coord2014, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2014$pop50km <- extract(pop_2014, coord2014, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2014$pop25km <- extract(pop_2014, coord2014, buffer = 25000, fun = sum)
toc()

# 2015

#sites with coords
#site name lat lon
coord2015 <- df2015[,c(3,5,6)]
head(coord2015)

#define coordiantes columns
coordinates(coord2015)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2015) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2015$pop100km <- extract(pop_2015, coord2015, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2015$pop50km <- extract(pop_2015, coord2015, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2015$pop25km <- extract(pop_2015, coord2015, buffer = 25000, fun = sum)
toc()

# 2016

#sites with coords
#site name lat lon
coord2016 <- df2016[,c(3,5,6)]
head(coord2016)

#define coordiantes columns
coordinates(coord2016)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2016) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2016$pop100km <- extract(pop_2016, coord2016, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2016$pop50km <- extract(pop_2016, coord2016, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2016$pop25km <- extract(pop_2016, coord2016, buffer = 25000, fun = sum)
toc()

# 2017

#sites with coords
#site name lat lon
coord2017 <- df2017[,c(3,5,6)]
head(coord2017)

#define coordiantes columns
coordinates(coord2017)<- c("Longitude.Degrees", "Latitude.Degrees")

#define spatial projections for sites
proj4string(coord2017) <- crs.geo


tic()
# Extract human population density for the world within 100 km buffer
df2017$pop100km <- extract(pop_2017, coord2017, buffer = 100000, fun = sum)

# Extract human population density for the world within 50 km buffer 
df2017$pop50km <- extract(pop_2017, coord2017, buffer = 50000, fun = sum)

# Extract human population density for the world within 25 km buffer 
df2017$pop25km <- extract(pop_2017, coord2017, buffer = 25000, fun = sum)
toc()


