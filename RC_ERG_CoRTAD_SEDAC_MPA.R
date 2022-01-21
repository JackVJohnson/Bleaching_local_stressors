rm(list=ls())

# overlay mpa data with bleaching and human pop

library(rgdal)

# mpa directory

mpashpwd<-"C:/Users/Jack V Johnson/Desktop/mpa shp"
setwd(mpashpwd)

mpafile <- readOGR(file.path(mpashpwd,"WDPA_Apr2020_marine-shapefile-polygons.shp"))

# bleaching directory

Bleaching_data_directory="C:/Users/Jack V Johnson/Desktop/Bleaching_synergistic"
setwd(Bleaching_data_directory)
Bleaching_Data <- read.csv("RC_ERG_SST_DHW_Pop.csv")


coordinates(Bleaching_Data)<- ~Longitude.Degrees+Latitude.Degrees
proj4string(Bleaching_Data)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
bldata<-spTransform(Bleaching_Data,proj4string(mpafile))

test<-over(bldata,mpafile)
head(test)

Bleaching_Data$MPA <- as.character(test$PA_DEF)
Bleaching_Data$STATUS_YR <- as.character(test$STATUS_YR)

Bleaching_Data <- as.data.frame(Bleaching_Data)
str(Bleaching_Data$MPA)
str(Bleaching_Data$STATUS_YR)

Bleaching_Data$MPA[Bleaching_Data$MPA > 0] <- "MPA"
Bleaching_Data$MPA[is.na(Bleaching_Data$MPA)] <- "Non MPA"

str(Bleaching_Data$MPA)

write.csv(Bleaching_Data, file = "RC_ERG_SST_DHW_Pop_MPA.csv")
