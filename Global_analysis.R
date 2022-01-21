
library(tidyverse)
library(MCMCglmm)
library(dotwhisker)
library(patchwork)
library(broom.mixed)
library(bayesplot)

rm(list=ls())

# bleaching directory

Bleaching_data_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/codes and files"
setwd(Bleaching_data_directory)
Bleaching_Data <- read.csv("Bleaching_data_new_pop.csv")

output_directory="C:/Users/40274182/OneDrive - Queen's University Belfast/PhD/Chapters/Bleaching/Local stressors/Figures and tables"
setwd(output_directory)

# tidy data

Bleaching_Data$Average_bleaching <- Bleaching_Data$Count / 4

Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$SSTA_DHW),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$Ecoregion),]

Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop25km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop50km),]
Bleaching_Data <- Bleaching_Data[!is.na(Bleaching_Data$pop100km),]

Bleaching_Data$human25km <- round(Bleaching_Data$pop25km, digits = 0)
Bleaching_Data$human50km <- round(Bleaching_Data$pop50km, digits = 0)
Bleaching_Data$human100km <- round(Bleaching_Data$pop100km, digits = 0)
Bleaching_Data$Average_bleaching <- round(Bleaching_Data$Average_bleaching, digits = 0)

summary(Bleaching_Data$human50km)
hist(Bleaching_Data$human50km)

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


Bleaching_Data$Ecoregion <- as.factor(Bleaching_Data$Ecoregion)
dev.off()
hist(Bleaching_Data$SSTA_DHW)
hist(Bleaching_Data$Average_bleaching)
hist(Bleaching_Data$human50km)

#  standardising 

model_df <- select(Bleaching_Data, SSTA_DHW, pop25km, pop50km, pop100km)

standardize_function<-function(x){
  x.standardized=(x-mean(na.omit(x)))/sd(na.omit(x))
  return(x.standardized)
}

model_df <- model_df; for(i in 1:ncol(model_df)) model_df[,i] <- standardize_function(model_df[,i])

temp <- select(Bleaching_Data, Average_bleaching, Ecoregion, STATUS_YR)

model_df <- cbind(model_df, temp)

hist(model_df$SSTA_DHW)
hist(model_df$pop25km)
hist(model_df$pop50km)
hist(model_df$pop100km)

# preliminary analysis 

mod_std_25 <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop25km, random = ~ Ecoregion, data = model_df, family = "poisson")
summary(mod_std_25)
plot(mod_std_25$Sol)
plot(mod_std_25$VCV)

mod_std_50 <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion, data = model_df, family = "poisson")
summary(mod_std_50)
plot(mod_std_50$Sol)
plot(mod_std_50$VCV)
hist(mcmc(mod_std_50$VCV)[,"Ecoregion"])

mod_std_100 <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop100km, random = ~ Ecoregion, data = model_df, family = "poisson")
summary(mod_std_100)

dic25 <-mod_std_25$DIC
dic50 <-mod_std_50$DIC
dic100 <-mod_std_100$DIC
DICs <- cbind(dic25, dic50, dic100)

output25 <- summary.MCMCglmm(mod_std_25)
output25 <- output25$solutions
output25 <- output25[-1,]

output50 <- summary.MCMCglmm(mod_std_50)
output50 <- output50$solutions
output50 <- output50[-1,]

output100 <- summary.MCMCglmm(mod_std_100)
output100 <- output100$solutions
output100 <- output100[-1,]
prelim <- rbind(output25, output50, output100)
prelim <- as.data.frame(prelim)

write.csv(DICs, file = "Model selection DICs.csv")
write.csv(prelim, file = "Preliminary coeffs.csv")


# formal analyses 

mod_std_50 <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion, data = model_df, family = "poisson")
summary(mod_std_50)
plot(mod_std_50$Sol)
plot(mod_std_50$VCV)
hist(mcmc(mod_std_50$VCV)[,"Ecoregion"])

# need to break the models up into mpa vs non mpa
# re-run model with paramater expanded prior for ecoregion 


global_model <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion, data = model_df, family = "poisson", nitt = 60000, burnin = 30000)
summary(global_model)
plot(global_model$Sol)
plot(global_model$VCV)
hist(mcmc(global_model$VCV)[,"Ecoregion"])

glob_output <- summary.MCMCglmm(global_model)
glob_output <- glob_output$solutions

mpa_sites <- subset(model_df, STATUS_YR == "MPA", select = 1:7)
non_mpa_sites <- subset(model_df, STATUS_YR == "Non MPA", select = 1:7)

mod_mpa <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion, data = mpa_sites, family = "poisson",  nitt = 60000, burnin = 30000)
summary(mod_mpa)
plot(mod_mpa$Sol)
plot(mod_mpa$VCV)

mpa_output <- summary.MCMCglmm(mod_mpa)
mpa_output <- mpa_output$solutions 

mod_nonmpa <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion, data = non_mpa_sites, family = "poisson", nitt = 60000, burnin = 30000)
summary(mod_nonmpa)
plot(mod_nonmpa$Sol)
plot(mod_nonmpa$VCV)

nonmpa_output <- summary.MCMCglmm(mod_nonmpa)
nonmpa_output <- nonmpa_output$solutions 

glob_output <- glob_output[-1,]
mpa_output <- mpa_output[-1,]
nonmpa_output <- nonmpa_output[-1,]

fig2_models <- rbind(glob_output, mpa_output, nonmpa_output)
fig2_models <- as.data.frame(fig2_models)

write.csv(fig2_models, file = "model_coeffs_fig2.csv")
library(dotwhisker)
library(tidybayes)
library(coda)
library(bayesplot)

dwplot(list(global_model, mod_mpa, mod_nonmpa))

tidyg <- broom.mixed::tidy(global_model) %>% 
  mutate(model = "global")
tidyg <- tidyg[-c(5,6),]
tidympa <- broom.mixed::tidy(mod_mpa)  %>% 
  mutate(model = "MPA")
tidympa <- tidympa[-c(5,6),]
tidynon <- broom.mixed::tidy(mod_nonmpa)  %>% 
  mutate(model = "non-MPA")
tidynon <- tidynon[-c(5,6),]

tidymods <- bind_rows(tidyg, tidympa, tidynon)
dwplot(tidymods)


tiff(file=file.path(output_directory,'model_coeffs_fig2.tif'),height=2000,width=3000,res=550)
dwplot(tidymods,
       vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>%
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  scale_color_viridis_d(name = "Model", 
                        labels = c("Global", "MPA", "Non-MPA"))
dev.off()

# trace plots 

tiff(file=file.path(output_directory,'Global_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(global_model$Sol)
dev.off()

tiff(file=file.path(output_directory,'MPA_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(mod_mpa$Sol)
dev.off()

tiff(file=file.path(output_directory,'NonMPA_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(mod_nonmpa$Sol)
dev.off()

tiff(file=file.path(output_directory,'Global_Eco_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(global_model$VCV)
dev.off()

tiff(file=file.path(output_directory,'MPA_Eco_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(mod_mpa$VCV)
dev.off()

tiff(file=file.path(output_directory,'nonMPA_Eco_Trace_plot.tif'),height=2000,width=2700,res=300)
plot(mod_nonmpa$VCV)
dev.off()

# ecoregion posterior variance

# global_model, mod_mpa, mod_nonmpa

mglob <- reshape2::melt(as.matrix(global_model$VCV[,"Ecoregion"]))
globhist <- ggplot(mglob,aes(value)) +
  geom_histogram(colour="black", binwidth = 0.5) + 
  scale_x_continuous(limits=c(2,11), breaks=seq(2,10, by=2)) +
  xlab("Posterior variance") + ylab("Frequency") +
  ggtitle("Global model") +
  theme_classic()

mmpa <- reshape2::melt(as.matrix(mod_mpa$VCV[,"Ecoregion"]))
mpahist <- ggplot(mmpa,aes(value))+
  ggtitle("MPA model") +
  scale_x_continuous(limits=c(2,11), breaks=seq(2,10, by=2)) +
  xlab("Posterior variance") + ylab("Frequency") +
  geom_histogram(colour="black", binwidth = 0.5) +
  theme_classic()

mnon <- reshape2::melt(as.matrix(mod_nonmpa$VCV[,"Ecoregion"]))
nonhist <- ggplot(mnon,aes(value))+
  ggtitle("Non-MPA model") +
  scale_x_continuous(limits=c(2,11), breaks=seq(2,10, by=2)) +
  xlab("Posterior variance") + ylab("Frequency") +
  geom_histogram(colour="black", binwidth = 0.5) +
  theme_classic()

tiff(file=file.path(output_directory,'Ecoregion_variance_fig3.tif'),height=2000,width=4000,res=400)
globhist + mpahist + nonhist + plot_layout(ncol = 3)
dev.off()


# mpa as a random effectt

str(model_df$STATUS_YR)
model_df$STATUS_YR <- as.factor(model_df$STATUS_YR)

mparandom_model <- MCMCglmm(fixed = Average_bleaching~SSTA_DHW*pop50km, random = ~ Ecoregion + STATUS_YR, data = model_df, family = "poisson", nitt = 60000, burnin = 30000)
summary(mparandom_model)
plot(mparandom_model$Sol)
plot(mparandom_model$VCV)

random_output <- tidy(mparandom_model)
random_output <- random_output[-c(5:7),]

tiff(file=file.path(output_directory,'mpa_as_random.tif'),height=2000,width=3000,res=550)
dwplot(random_output, 
       vline = geom_vline(xintercept = 0, colour = "grey", linetype = 2)) %>%
  relabel_predictors(c(SSTA_DHW = "DHW",
                       pop50km = "Local Population",
                       "SSTA_DHW:pop50km" = "Interaction")) + # as character string bevause of stupid : 
  theme_classic() + 
  xlab("Coefficient Estimate") + 
  ylab("") +
  scale_color_viridis_d() +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
dev.off()
