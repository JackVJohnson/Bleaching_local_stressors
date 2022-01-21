
library(tidyverse)
library(rgdal)
library(RColorBrewer)
library(plotrix)
#install.packages('GISTools')
library(GISTools)
library(R2OpenBUGS)
#install.packages('viridis')
library(viridis)
#install.packages('Cairo')
library(Cairo)

#working directory 

worldmapwd<-"C:/Users/Jack V Johnson/Desktop/World map"
setwd(worldmapwd)
wlrd.p <- readOGR(file.path(worldmapwd,'TM_WORLD_BORDERS_SIMPL_PC150.shp'))

Bleaching_data_directory="C:/Users/Jack V Johnson/Desktop/Bleaching_synergistic"
Output_directory="C:/Users/Jack V Johnson/Desktop/Bleaching_synergistic/Output"

setwd(Bleaching_data_directory)
Bleaching_Data <- read.csv("RC_ERG_SST_DHW_Pop_MPA.csv")

Bleaching_Data$Average_bleaching <- Bleaching_Data$Count / 4

Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$SSTA_DHW),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$Ecoregion),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$human25km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$human50km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$human100km),]

Bleaching_Data$human25km <- round(Bleaching_Data$human25km, digits = 0)
Bleaching_Data$human50km <- round(Bleaching_Data$human50km, digits = 0)
Bleaching_Data$human100km <- round(Bleaching_Data$human100km, digits = 0)
Bleaching_Data$Average_bleaching <- round(Bleaching_Data$Average_bleaching, digits = 0)

Bleaching_Data$STATUS_YR[Bleaching_Data$STATUS_YR>Bleaching_Data$Year] <- NA
Bleaching_Data$STATUS_YR[Bleaching_Data$STATUS_YR>1] <- "MPA"
Bleaching_Data$STATUS_YR[is.na(Bleaching_Data$STATUS_YR)] <- "Non MPA"
Bleaching_Data$STATUS_YR

Bleachmpa=subset(Bleaching_Data, STATUS_YR == "MPA")
Bleachnon=subset(Bleaching_Data, STATUS_YR == "Non MPA")

Bleachnon$freq <- 1
dfnon <- aggregate(Bleachnon$freq, by=list(Site=Bleachnon$Site, Lat=Bleachnon$Latitude.Degrees, Lon=Bleachnon$Longitude.Degrees), FUN=sum)

Bleachmpa$freq <- 1
dfmpa <- aggregate(Bleachmpa$freq, by=list(Site=Bleachmpa$Site, Lat=Bleachmpa$Latitude.Degrees, Lon=Bleachmpa$Longitude.Degrees), FUN=sum)


# make the map

windowsFonts(Arial=windowsFont("TT Arial"))

par(family="Arial")
pal<-viridis(7)
tiff(file=file.path(Output_directory,'nonmp_-mpa_surveys.tiff'),height=1800,width=3500,res=300)
par(mfrow=c(2,1), mgp=c(0.5,0.6,0), mar=c(1,1,1,1))
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(-2000000,2000000), col='grey90',border='grey70')
title(main="(a)",adj=0)
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
box()
xy <- dfnon[dfnon$x == 0,c('Lon','Lat')]
xy <- SpatialPointsDataFrame(dfnon=xy,coords=xy[c('Lon','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)
temp <- subset(dfnon, x > 0)
temp <- temp[with(temp, order(temp$x)),]
xy <- temp[c('Lon','Lat')]
xy <- SpatialPointsDataFrame(data=xy,coords=xy[c('Lon','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)
points(xy, cex=.7, pch=19, col=pal[temp$x])
text(-7868896,-2922012,'Indian Ocean',cex=.8)
text(9438742,487176,'Pacific Ocean',cex=.8)
north.arrow(x=(-16654136+111319.4*320), y=1615153*2, len=(111319.4*2), lab="N", cex=.7)
#legend
plotrix::color.legend(9684797.171+25e5,-28*111319.4666666667,15807371.62+25e5,-23.5*111319.4666666667,legend=c(0,43),rect.col=c("white",pal),cex=1)
text(((15807371.62+25e5)-(9684797.171+25e5))/2+(9684797.171+25e5),-18*111319.4666666667,"Number of Surveys", cex=.75)

# mpa
plot(wlrd.p,ylim=c(-4400000,4400000),xlim=c(-2000000,2000000), col='grey90',border='grey70')
title(main="(b)",adj=0)
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
box()
xy <- dfmpa[dfmpa$x == 0,c('Lon','Lat')]
xy <- SpatialPointsDataFrame(dfman=xy,coords=xy[c('Lon','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)
temp <- subset(dfmpa, x > 0)
temp <- temp[with(temp, order(temp$x)),]
xy <- temp[c('Lon','Lat')]
xy <- SpatialPointsDataFrame(data=xy,coords=xy[c('Lon','Lat')], proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
xy <- spTransform(xy,CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=150 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points(xy, cex=.7)
points(xy, cex=.7, pch=19, col=pal[temp$x])
text(-7868896,-2922012,'Indian Ocean',cex=.8)
text(9438742,487176,'Pacific Ocean',cex=.8)
north.arrow(x=(-16654136+111319.4*320), y=1615153*2, len=(111319.4*2), lab="N", cex=.7)
#legend
plotrix::color.legend(9684797.171+25e5,-28*111319.4666666667,15807371.62+25e5,-23.5*111319.4666666667,legend=c(0,43),rect.col=c("white",pal),cex=1)
text(((15807371.62+25e5)-(9684797.171+25e5))/2+(9684797.171+25e5),-18*111319.4666666667,"Number of Surveys", cex=.75)

dev.off()