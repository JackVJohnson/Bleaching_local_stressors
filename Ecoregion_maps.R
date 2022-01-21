library(MASS)
library(audio)
library(sp)
library(foreign)
library(rgdal)
library(maptools)
library(rgeos)
library(doParallel)
library(rasterVis)
library(dismo)
library(plotKML)
library(PBSmapping)
library(lme4)
library(blme)
library(raster)
library(RColorBrewer)
library(sjmisc)
library(ncdf4)
library(knitr)
library(stringr)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(dplyr)
library(maptools)
library(plotrix)
library(GISTools)
library(viridis)
library(fields)
library(dichromat)
library(colorspace)
library(maps)

Bleaching_data_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/codes and files"
Output_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/Figures and tables"
shp_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/ERG_shp_files"


ecoregions <- readOGR(file.path(shp_directory,'ecoregion_exportPolygon.shp'))


ecos_list<-c()
for (i in 1:150){
  eco_i<-Polygons((Filter(function(f){f@ringDir==1}, ecoregions@polygons[[i]]@Polygons)), ID=i)
  ecos_list<-append(ecos_list, values=eco_i, after = length(ecos_list))
  Sys.sleep(.2)
  print(i)
}
ecos<-SpatialPolygons(ecos_list)

ecos$ERG<-ecoregions$ERG
ecos$Ecoregion<-ecoregions$Ecoregion
ecos@proj4string<-ecoregions@proj4string
ecos@plotOrder<-ecoregions@plotOrder
ecos@data<-ecoregions@data
ecoregions<-ecos



plot.map<- function(database,center,transf=T,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  newproj <- "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #utm
  nextproj<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" #latlong
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  
  
  colnames(polygons)<-c("x",'y')
  polygons<-as.data.frame(polygons)
  z<-complete.cases(polygons)
  p<-z
  z<-cbind(z,z)
  polygons<-polygons[complete.cases(polygons),]
  coordinates(polygons)<-~x+y
  proj4string(polygons)<-CRS(nextproj)
  if(transf==T){ polygons<-spTransform(polygons,CRS(newproj))}
  
  z[p==F,]<-c(NA,NA)
  z[which(p==T),]<-coordinates(polygons)
  Obj[[1]] <- z[,1]
  Obj[[2]] <- z[,2]
  
  map(Obj,...)
}


setwd(Bleaching_data_directory)
getwd()
ecoregions_probs<-read.csv(file = "ecoregion_p_values.csv", header = TRUE, sep=",")

names(ecoregions_probs)
colnames(ecoregions_probs)[1] <- "ERG"
ecoregions_probs$ERG <- as.character(ecoregions_probs$ERG)
#cols <- RColorBrewer::brewer.pal(n=5, "Spectral")
cols <- colorRampPalette(colorschemes$DarkRedtoBlue.12)(5)
#cols <- diverging_hcl(5, "Blue-Red 3")
brks<-seq(from=-3, to=2, by=1)

ecoregions_probs$glob_dhw_cols <- cols[(cut(ecoregions_probs$dhw_glob, brks))]
ecoregions_probs$glob_pop_cols <- cols[(cut(ecoregions_probs$population_glob, brks))]
ecoregions_probs$glob_int_cols <- cols[(cut(ecoregions_probs$interaction_glob, brks))]


ecoregions$glob_dhw_cols<-"white"
ecoregions$glob_pop_cols<-"white"
ecoregions$glob_int_cols<-"white"

for (i in 1:150){
  ecoregions$glob_dhw_cols[i]<-ecoregions_probs$glob_dhw_cols[ecoregions$ERG[i]==ecoregions_probs$ERG]
  ecoregions$glob_pop_cols[i]<-ecoregions_probs$glob_pop_cols[ecoregions$ERG[i]==ecoregions_probs$ERG]
  ecoregions$glob_int_cols[i]<-ecoregions_probs$glob_int_cols[ecoregions$ERG[i]==ecoregions_probs$ERG]
}

setwd(shp_directory)

tiff(file=file.path(Output_directory, 'Ecoregion_maps.tif'), ,height=3000,width=3800,res=300)
#par(mgp=c(0.5,0.6,0), mar=c(0,0,0,0))
par(mfrow=c(3,1), mgp=c(0.5,0.6,0), mar=c(0,0,0,0))

# global dhw 
plot(ecoregions, col=ecoregions$glob_dhw_cols, border="turquoise4", lwd=1, xlim=c(111319.4*-150,111319.4*150), ylim=c(111319.4*-43, 111319.4*50))
box()
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=1.0, family='A')
text(9438742,487176,'Pacific Ocean',cex=1.0, family='A')
text(x=(-16654136+111319.4*305), y=1615153*2.15,'Atlantic Ocean',cex=1.0, family='A')
north.arrow(x=(-16654136+111319.4*220), y=1615153*3, len=(111319.4*2), lab="N", cex=.7)
scalebar(d=111319.4*36, xy=c((-16654136+111319.4*-27), (1615153*-2.7)), type='bar', divs=4, below="kilometers", label=c("0", "", "4000"))
mapcol  <- cols
plotrix::color.legend(60*111319.4,-43*111319.4,120*111319.4,-38*111319.4,legend=c(""),rect.col=cols,cex=1)
text(60*111319.4,-35*111319.4,'-')
text(120*111319.4,-35*111319.4,'+')
text(90*111319.4,-35*111319.4,'Bleaching Probability')
plot.map("world", center=0 ,bg="#00000000",fill=TRUE, col="lightgray",add=T,xlab='longitude',ylab='latitude')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
text(-20654136, 5461227, '(a)', cex = 1.5, font =2)

# global pop

plot(ecoregions, col=ecoregions$glob_pop_cols, border="turquoise4", lwd=1, xlim=c(111319.4*-150,111319.4*150), ylim=c(111319.4*-43, 111319.4*50))
box()
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=1.0, family='A')
text(9438742,487176,'Pacific Ocean',cex=1.0, family='A')
text(x=(-16654136+111319.4*305), y=1615153*2.15,'Atlantic Ocean',cex=1.0, family='A')
north.arrow(x=(-16654136+111319.4*220), y=1615153*3, len=(111319.4*2), lab="N", cex=.7)
scalebar(d=111319.4*36, xy=c((-16654136+111319.4*-27), (1615153*-2.7)), type='bar', divs=4, below="kilometers", label=c("0", "", "4000"))
mapcol  <- cols
plotrix::color.legend(60*111319.4,-43*111319.4,120*111319.4,-38*111319.4,legend=c(""),rect.col=cols,cex=1)
text(60*111319.4,-35*111319.4,'-')
text(120*111319.4,-35*111319.4,'+')
text(90*111319.4,-35*111319.4,'Bleaching Probability')
plot.map("world", center=0 ,bg="#00000000",fill=TRUE, col="lightgray",add=T,xlab='longitude',ylab='latitude')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
text(-20654136, 5461227, '(b)', cex = 1.5, font =2)

# global interaction

plot(ecoregions, col=ecoregions$glob_int_cols, border="turquoise4", lwd=1, xlim=c(111319.4*-150,111319.4*150), ylim=c(111319.4*-43, 111319.4*50))
box()
windowsFonts(A=windowsFont('Arial Unicode MS'))
text(-7868896,-2922012,'Indian Ocean',cex=1.0, family='A')
text(9438742,487176,'Pacific Ocean',cex=1.0, family='A')
text(x=(-16654136+111319.4*305), y=1615153*2.15,'Atlantic Ocean',cex=1.0, family='A')
north.arrow(x=(-16654136+111319.4*220), y=1615153*3, len=(111319.4*2), lab="N", cex=.7)
scalebar(d=111319.4*36, xy=c((-16654136+111319.4*-27), (1615153*-2.7)), type='bar', divs=4, below="kilometers", label=c("0", "", "4000"))
mapcol  <- cols
plotrix::color.legend(60*111319.4,-43*111319.4,120*111319.4,-38*111319.4,legend=c(""),rect.col=cols,cex=1)
text(60*111319.4,-35*111319.4,'-')
text(120*111319.4,-35*111319.4,'+')
text(90*111319.4,-35*111319.4,'Bleaching Probability')
plot.map("world", center=0 ,bg="#00000000",fill=TRUE, col="lightgray",add=T,xlab='longitude',ylab='latitude')
axis(1,at=c(-10018754.17,3339584.724,16697920),lab=c('60°','180°','-60° '),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(2, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('23°','0°','-23°'),las=3,tcl=0.35,mgp=c(-2,-1.3,0),hadj=.4)
axis(3,at=c(-10018754.17,3339584.724,16697920),lab=c('','',''),las=1,tcl=0.35,mgp=c(-1,-1.3,0))
axis(4, at=c(23*111319.4666666667,0,-23*111319.4666666667),labels=c('','',''),las=2,tcl=0.35,mgp=c(-1,-0.6,0),hadj=0)
text(-20654136, 5461227, '(c)',cex = 1.5, font =2)

dev.off()