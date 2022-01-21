# analysis by ecoregion

library(rgdal)
library(MCMCglmm)
library(tidyverse)
library(dotwhisker)
library(broom.mixed)

Bleaching_data_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/codes and files"
setwd(Bleaching_data_directory)
Bleaching_Data <- read.csv("Bleaching_data_new_pop.csv")

output_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/Figures and tables"

Ecoregion_directory <- "C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/ERG_shp_files"
setwd(Ecoregion_directory)

# open ecoregions to add ecoregion names to data 

ECO<-readOGR("ecoregion_exportPolygon.shp")

# get rid of holes in shp 

ecos_list<-c()
for (i in 1:150){
  eco_i<-Polygons((Filter(function(f){f@ringDir==1}, ECO@polygons[[i]]@Polygons)), ID=i)
  ecos_list<-append(ecos_list, values=eco_i, after = length(ecos_list))
  #include a brief pause (the Sys.sleep function) because if running in Rstudio, it takes a while for the code to run and for the value to be loaded into the global environment. If there is no pause, the next iteration of the loop starts before the previous value is fully saved and loaded into the environment, and there can be errors in the shapefiles 
  Sys.sleep(.2)
}

ecos<-SpatialPolygons(ecos_list)
ecos$Ecoregion<-ECO$Ecoregion
ecos@proj4string<-ECO@proj4string
ecos@plotOrder<-ECO@plotOrder
ecos@data<-ECO@data
ECO<-ecos

# overlay 

coordinates(Bleaching_Data)<- ~Longitude.Degrees+Latitude.Degrees
proj4string(Bleaching_Data)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
bldata<-spTransform(Bleaching_Data,proj4string(ECO))

test<-over(bldata,ECO)
Bleaching_Data$Ecoregion_name <- as.character(test$Ecoregion)

Bleaching_Data <- as.data.frame(Bleaching_Data)

# tidy data 

Bleaching_Data$Average_bleaching <- Bleaching_Data$Count / 4

Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$SSTA_DHW),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$Ecoregion),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop25km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop50km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop100km),]

Bleaching_Data$pop25km <- round(Bleaching_Data$pop25km, digits = 0)
Bleaching_Data$pop50km <- round(Bleaching_Data$pop50km, digits = 0)
Bleaching_Data$pop100km <- round(Bleaching_Data$pop100km, digits = 0)
Bleaching_Data$Average_bleaching <- round(Bleaching_Data$Average_bleaching, digits = 0)

summary(Bleaching_Data$pop50km)
hist(Bleaching_Data$pop50km)

Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$Average_bleaching),]
hist(Bleaching_Data$Average_bleaching)


# mpa based on year of designation

summary(Bleaching_Data$STATUS_YR)
head(Bleaching_Data$STATUS_YR)
tail(Bleaching_Data$STATUS_YR)
#mtransform mpa Bleaching_Dataata which are NA to workable format specific for each column - NA signifies site Bleaching_Dataoes not fall within an mpa, therefore shoulBleaching_Data be analyseBleaching_Data as a non protecteBleaching_Data area

Bleaching_Data$STATUS_YR[Bleaching_Data$STATUS_YR>Bleaching_Data$Year] <- NA
Bleaching_Data$STATUS_YR[Bleaching_Data$STATUS_YR>1] <- "MPA"
Bleaching_Data$STATUS_YR[is.na(Bleaching_Data$STATUS_YR)] <- "Non MPA"
Bleaching_Data$STATUS_YR
##########

model_df <- select(Bleaching_Data, SSTA_DHW, pop25km, pop50km, pop100km)

standardize_function<-function(x){
  x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
  return(x.standardized)
}

model_df <- model_df; for(i in 1:ncol(model_df)) model_df[,i] <- standardize_function(model_df[,i])

temp <- select(Bleaching_Data, Average_bleaching, Ecoregion, STATUS_YR, Ecoregion_name)

model_df <- cbind(model_df, temp)

hist(model_df$pop50km)

########

model_df_mpa <- subset(Bleaching_Data, STATUS_YR == "MPA", select = 1:8)
model_df_nonmpa <- subset(Bleaching_Data, STATUS_YR == "Non MPA", select = 1:8)

freqofeco <- data.frame(Bleaching_Data$Ecoregion_name)
freqofeco <- as.data.frame(table(freqofeco))
eco_30plus <- subset(freqofeco, Freq > 30)
print(eco_30plus$freqofeco)


# run models for each eco region
options(scipen = 999)

modAND <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Andaman Sea'), family = "poisson", nitt = 1500000, burnin = 900000)
summary(modAND) 

modAND_out <- tidy(modAND) %>%
  filter(!row_number()==5)

modBFK <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Bahamas and Florida Keys'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBFK)

modBFK_out <- tidy(modBFK) %>%
  filter(!row_number()==5)

modBMI <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Banda Sea and Molucca Islands'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modBMI) # discard
plot(modBMI$Sol)

modBWC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Belize and west Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBWC)

modBWC_out <- tidy(modBWC) %>%
  filter(!row_number()==5)

modBNG <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Bismarck Sea, New Guinea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modBNG) # discard

modBRA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Brazil'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBRA)

modBRA_out <- tidy(modBRA) %>%
  filter(!row_number()==5)

modCEL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Celebes Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modCEL)

modCEL_out <- tidy(modCEL) %>%
  filter(!row_number()==5)

modGBR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Central and northern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGBR)

modGBR_out <- tidy(modGBR) %>%
  filter(!row_number()==5)

modCKA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Cocos Keeling Atolls, Indian Ocean'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modCKA) # discard

modEHA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Eastern Hawaii'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modEHA)

modEHA_out <- tidy(modEHA) %>%
  filter(!row_number()==5)

modFIJ <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Fiji'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modFIJ)

modFIJ_out <- tidy(modFIJ) %>%
  filter(!row_number()==5)

modGOO <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Gulf of Oman'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGOO)
plot(modGOO$Sol)

modGOO_out <- tidy(modGOO) %>%
  filter(!row_number()==5)

modGOT <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Gulf of Thailand'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modGOT)
plot(modGOT$Sol) 

modGOT_out <- tidy(modGOT) %>%
  filter(!row_number()==5)

modHPR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Hispaniola, Puerto Rico and Lesser Antilles'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHPR)

modHPR_out <- tidy(modHPR) %>%
  filter(!row_number()==5)

modHKN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Hong Kong'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHKN)

modHKN_out <- tidy(modHKN) %>%
  filter(!row_number()==5)

modJAM <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Jamaica'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modJAM)

modJAM_out <- tidy(modJAM) %>%
  filter(!row_number()==5)

modJAV <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Java Sea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modJAV) # discard

modLSI <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Lesser Sunda Islands and Savu Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modLSI)

modLSI_out <- tidy(modLSI) %>%
  filter(!row_number()==5)

modMSI <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Makassar Strait, Indonesia'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modMSI) # discard
plot(modMSI$Sol)

modMAL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Maldive Islands'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMAL)

modMAL_out <- tidy(modMAL) %>%
  filter(!row_number()==5)

modMSC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Mascarene Islands'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modMSC) # discard

modMOC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Mayotte and Comoros'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMOC) # despite convergence, values in 1000s so discard (also not significant so highly misleading) 
plot(modMOC$Sol)

modMBE <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Moreton Bay, eastern Australia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMBE)

modMBE_out <- tidy(modMBE) %>%
  filter(!row_number()==5)

modNES <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Netherlands Antilles and south Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNES)

modNES_out <- tidy(modNES) %>%
  filter(!row_number()==5)

modNCD <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='New Caledonia'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modNCD) # discard

modNCR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='North and central Red Sea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modNCR) # discard

modNPH <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='North Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNPH) 

modNPH_out <- tidy(modNPH) %>%
  filter(!row_number()==5)

modVIE <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='North Vietnam'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modVIE) # discard

modPKM <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Pohnpei and Kosrae, Micronesia'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modPKM) # discard

modSIF <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Society Islands, French Polynesia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSIF)

modSIF_out <- tidy(modSIF) %>%
  filter(!row_number()==5)

modSIB <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Solomon Islands and Bougainville'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSIB) # discard

modSEP <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='South-east Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSEP)

modSEP_out <- tidy(modSEP) %>%
  filter(!row_number()==5)

modSRJ <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='South Ryukyu Islands, Japan'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSRJ) # discard

modSVT <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='South Vietnam'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSVT) # discard

modSGB <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Southern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSGB)

modSGB_out <- tidy(modSGB) %>%
  filter(!row_number()==5)

modSUS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Sulu Sea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSUS)

modSUS_out <- tidy(modSUS) %>%
  filter(!row_number()==5)

modSSS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Sunda Shelf, south-east Asia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSSS)

modSSS_out <- tidy(modSSS) %>%
  filter(!row_number()==5)

modTCC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Taiwan and coastal China'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modTCC)

modTCC_out <- tidy(modTCC) %>%
  filter(!row_number()==5)

modVAN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Vanuatu'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modVAN)

modVAN_out <- tidy(modVAN) %>%
  filter(!row_number()==5)

modWTA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df,Ecoregion_name=='Western Tuamotu Archipelago, central Pacific'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modWTA) # discard 

# ECOREGION COEFFICIENT PLOTS


P1 <- dwplot(modAND_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Andaman Sea") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


P2 <- dwplot(modBFK_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  scale_color_viridis_d() +
  ggtitle("Bahamas and Florida Keys") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


P3 <- dwplot(modBWC_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  scale_color_viridis_d() +
  ggtitle("Belize and West Caribbean") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


P4 <- dwplot(modBRA_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Brazil") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P5 <- dwplot(modCEL_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Celebes Sea") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P6<- dwplot(modEHA_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Eastern Hawaii") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P7 <- dwplot(modFIJ_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Fiji") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P8 <- dwplot(modGBR_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Central and Northern GBR") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P9 <- dwplot(modGOO_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Gulf of Oman") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P10 <- dwplot(modGOT_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Gulf of Thailand") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P11 <- dwplot(modHKN_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Hong Kong") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P12 <- dwplot(modHPR_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Hispaniola") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P13 <- dwplot(modJAM_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Jamaica") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P14 <- dwplot(modLSI_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Lesser Sunda Islands and Savi Sea") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P15 <- dwplot(modMAL_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Maldives") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P16 <- dwplot(modMBE_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Moreton Bay, Eastern Australia") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P17 <- dwplot(modNES_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Netherlands Antilles") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P18 <- dwplot(modNPH_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("North Philippines")+
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P19 <- dwplot(modSEP_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("South-East Philippines") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P20 <- dwplot(modSGB_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Southern GBR") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


P21 <- dwplot(modSIF_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Society Islands, French Polynesia") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P22 <- dwplot(modSSS_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Sunda Shelf, south-east Asia") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P23 <- dwplot(modSUS_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Sulu Sea") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P24 <- dwplot(modTCC_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Taiwan and coastal china") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")

P25 <- dwplot(modVAN_out, vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>% 
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  ggtitle("Vanuatu") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


library(patchwork)
setwd(output_directory)

png(file='ecoregion_coeff_plots.png',height=5000,width=6000,res=300)
P1+P2+P3+P4+P5+P6+P7+P8+P9+P10+P11+P12+P13+P14+P15+P16+P17+P18+P19+P20+P21+P22+P23+P24+P25+
  plot_layout(ncol = 5, nrow = 5)
dev.off()

#png(file='ecoregion_coeff_plots.png',height=5000,width=5000,res=300)
#(P1|P2|P3|P4)/(P5|P6|P7|P8)/(P9|P10|P11|P12)/(P13|P14|P15|P16)/(P17|P18|P19|P20)/(P21|P22|P23|P24)
#dev.off()




#####################################################

# mpa models

modAND <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Andaman Sea'), family = "poisson", nitt = 1500000, burnin = 900000)
summary(modAND) 

modBFK <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Bahamas and Florida Keys'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBFK)

modBWC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Belize and west Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBWC)

modBRA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Brazil'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBRA)

modCEL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Celebes Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modCEL)

modGBR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Central and northern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGBR)

modEHA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Eastern Hawaii'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modEHA)

modFIJ <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Fiji'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modFIJ)

modGOO <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Gulf of Oman'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGOO)
plot(modGOO$Sol)

modGOT <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Gulf of Thailand'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modGOT)
plot(modGOT$Sol) 

modHPR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Hispaniola, Puerto Rico and Lesser Antilles'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHPR)

modHKN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Hong Kong'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHKN)  

modJAM <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Jamaica'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modJAM)

modLSI <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Lesser Sunda Islands and Savu Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modLSI)

modMAL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Maldive Islands'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMAL) # discard

modMOC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Mayotte and Comoros'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMOC) # discard

modMBE <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Moreton Bay, eastern Australia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMBE)

modNES <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Netherlands Antilles and south Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNES)

modNPH <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='North Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNPH) # discard 

modSIF <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Society Islands, French Polynesia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSIF) # discard

modSEP <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='South-east Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSEP)

modSGB <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Southern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSGB)

modSUS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Sulu Sea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSUS)

modSSS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Sunda Shelf, south-east Asia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSSS)

modTCC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Taiwan and coastal China'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modTCC) # discard

modVAN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_mpa,Ecoregion_name=='Vanuatu'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modVAN) #discard

#########################

# non mpa sitees


modAND <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Andaman Sea'), family = "poisson", nitt = 1500000, burnin = 900000)
summary(modAND) # discard

modBFK <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Bahamas and Florida Keys'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBFK)

modBWC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Belize and west Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBWC)

modBRA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Brazil'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modBRA)

modCEL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Celebes Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modCEL)

modGBR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Central and northern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGBR) # discard 

modEHA <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Eastern Hawaii'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modEHA) 

modFIJ <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Fiji'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modFIJ)

modGOO <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Gulf of Oman'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modGOO)
plot(modGOO$Sol)

modGOT <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Gulf of Thailand'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modGOT)
plot(modGOT$Sol) 

modHPR <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Hispaniola, Puerto Rico and Lesser Antilles'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHPR)

modHKN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Hong Kong'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modHKN)  

modJAM <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Jamaica'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modJAM)

modLSI <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Lesser Sunda Islands and Savu Sea'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modLSI)

modMAL <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Maldive Islands'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMAL) 

modMOC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Mayotte and Comoros'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMOC) # discard

modMBE <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Moreton Bay, eastern Australia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modMBE)

modNES <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Netherlands Antilles and south Caribbean'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNES)

modNPH <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='North Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modNPH)

modSIF <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Society Islands, French Polynesia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSIF)

modSEP <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='South-east Philippines'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSEP)

modSGB <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Southern Great Barrier Reef'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSGB) # discard 

modSUS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Sulu Sea'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modSUS)

modSSS <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Sunda Shelf, south-east Asia'), family = "poisson", nitt = 100000, burnin = 50000)
summary(modSSS)

modTCC <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Taiwan and coastal China'), family = "poisson",nitt = 1000000, burnin = 600000)
summary(modTCC) 

modVAN <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, data = subset(model_df_nonmpa,Ecoregion_name=='Vanuatu'), family = "poisson", nitt = 1000000, burnin = 600000)
summary(modVAN)



