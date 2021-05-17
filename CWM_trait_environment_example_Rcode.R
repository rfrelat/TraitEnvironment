############################################################
# FutureMares Task 1.2. Trait-environment workshop. 3rd June 2021
# Example data and R code for estimating CWM traits and diversity indices, 
# as well as plotting patterns and identifying key drivers and state-pressure relationships (using GAMs and random forest)
# Martin Lindegren, DTU Aqua, 12/5 2021. (Includes also modifed code from Esther Beukhof - thanks!)
##############################################################
# Clear memory
rm(list=ls())

## Install packages
install.packages("FD")
install.packages("mgcv")
install.packages("qqplot2")
install.packages("RColorBrewer")
install.packages("randomForest")
install.packages("MuMIn")

# Load packages
library(FD) ## Calculate functional diversity indices
library(mgcv)
library(ggplot2)
library(RColorBrewer)
library(randomForest)
library(MuMIn)

######################################
# 1. Load and inspect data
##############################3
# Load data
setwd("H:/Work/Project coordination/FutureMares/WP coordination/Trait-env workshop")
load(file='NorthSea_FishTraitEnv.RData') # 

# Inspect data tables 
head(abu) # Species abundance by grid cell (the community matrix)
head(env) # Environmental conditions by grid cell
head(trait) # Traits by species
head(coo) # The spatial coordinates (lat/lon) for each grid cell

#########################################################
# 2. Calculate Community-weighted mean (CWM) traits and various diversity indices 
#########################################################
?dbFD  ## Inspect the function and all its option
FD<-dbFD(trait, abu, stand.FRic=T)  #stand.FRic=T, means FRic will range between 0 and 1

# Inspect results
FD  #### All results
FD$nbsp  ### Species Richness              
FD$FRic  ### Functional Richness           
FD$FEve  ### Functional Eveness            
FD$FDiv  ### Functional Divergence         
FD$RaoQ  ### Functional Entropy            
FD$CWM   ### Community weighted mean traits 

# Reformat into an output data frame
FD<-do.call(cbind, FD)

## Add lon, lat to output data frame
FD[,c("lon","lat")]<-coo

########################################
# 3. Plot biodiversity indices by lat,long
####################################
colpal<-rev(brewer.pal(11,"RdYlBu")) # Color palette for the map
# Plotting function 
funPlot<-function(Var){
  colnames(FD)[which(colnames(FD)==Var)]<-"Var"
  ggplot() +
    geom_tile(data=FD,aes(x=lon,y=lat,fill=Var)) + # 
    scale_fill_gradientn(name = Var, colours=colpal, na.value = 'white')+#,limits=c(0,25100)) +
    borders(fill="gray44",colour="black") +
    coord_quickmap(xlim=c(range(FD$lon)),ylim=c(range(FD$lat)))+
    labs(x = "Lon",y="Lat")+
    theme(legend.position = c(0.93,0.8))+
    theme(legend.title = element_text(size = 8),legend.text = element_text(size = 7))+
    guides(shape = guide_legend(override.aes = list(size = 0.2)))+
    theme(panel.background = element_rect(fill=alpha('light blue', 0.4), colour = 'black'))
}
Metrics<-colnames(FD)
funPlot(Var=Metrics[9])# Plot CWM traits and diversity metrics


#################################################################
# 4. Statistically investigate trait-environment relationships (i.e., here trait responses at the community level)
##############################################################
# Merge FD with env data frame
FDenv<-cbind(FD,env)
head(FDenv)
  
# Initial exploration of model fit for selected traits using all covariates 
# Test GAM
?gam # See documentation

funGam<-function(Var){
  colnames(FDenv)[which(colnames(FDenv)==Var)]<-"Var"
  fit<-gam(Var~s(Depth, k=3)+s(SBT,k=3)+s(SBS,k=3)+s(Chl,k=3)+s(SBT_sea,k=3)+ s(Chl_sea, k=3)+s(Fishing, k=3), na.action=na.exclude, data=FDenv)
  par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
  plot(fit,shade=T,shade.col="grey",res=T,rug=F,pch=20,ylab="")
  par(mfrow=c(1,1),mar=c(5,4,2,2))
  return(summary(fit))
}
#Plot GAM partial smooth plots and return summary statistics for the fitted model
funGam(Metrics[9])

#############################################
# More detailed and thorough model selection 
# Here we use a mixed GAM approach taking into account potential spatial autocorrelation (by allowing for correlated error structures)
?uGamm # Wrapper function for gamm - see specification

# Create dummy variable to include as random factor
dummy <- rep(1, dim(FDenv)[1]) ### Dummy variables to be add as random factor (mandatory, but won't 'change' anything)

# Trophic level as an example
# Check distribution of response
hist(FDenv$CWM.Trophic.level) # Consider log transforming or change family statement in gam in case of strong deviations from normal distribution 

fitTL <- uGamm(CWM.Trophic.level ~ s(Depth, k=3)
                   +s(SBT,k=3)
                   +s(SBS,k=3)
                   +s(Chl,k=3)
                   +s(SBT_sea,k=3)
                   +s(Chl_sea, k=3)
                   +s(Fishing, k=3),
                   gaussian(link = "identity"),random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat), # Can try other error structures
                   data = FDenv, control=lmeControl(opt="optim"))
# Inspect summary statistics
summary(fitTL$lme)
summary(fitTL$gam)

# Check diagnostics
par(mfrow=c(2,2))
gam.check(fitTL$gam)
par(mfrow=c(1,1))

# Plot model smooth terms
par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
plot(fitTL$gam,shade=T,shade.col="grey",res=T,rug=F,pch=20)
par(mfrow=c(1,1),mar=c(5,4,2,2))

# Perform automated model selection testing all combinations of predictors using the dredge function (and rank according to AICc)
?dredge
results <- dredge(fitTL, m.lim=c(1,4), rank="AICc") # Here test max 4 variables per model to reduce run time
subset(results, delta <5)  # Depth, SBT and fishing are key variables
# Calculate and view relative variable importance (RVI) scores
importance(results) 

# Fit and inspect the "best" model
fitTLb <- uGamm(CWM.Trophic.level ~ s(Depth, k=3)
               +s(SBT,k=3),
               #+s(Fishing, k=3),
               gaussian(link = "identity"),random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat), # Can try other error structures
               data = FDenv, control=lmeControl(opt="optim"))
# Inspect summary statistics
summary(fitTLb$lme)
summary(fitTLb$gam)

# Check diagnostics
par(mfrow=c(2,2))
gam.check(fitTLb$gam)
par(mfrow=c(1,1))

# Plot model smooth terms
par(mfrow=c(2,1),mar=c(4,4,0.2,0.2))
plot(fitTLb$gam,shade=T,shade.col="grey",res=T,rug=F,pch=20)
par(mfrow=c(1,1),mar=c(5,4,2,2))

###################################
# Test Random forest
?randomForest # See help file
 
# Wrapper function to explore RF for a given trait
funRf<-function(Var){
  colnames(FDenv)[which(colnames(FDenv)==Var)]<-"Var"
  fit<-randomForest(Var~Depth+SBT+SBS+Chl+SBT_sea+Chl_sea+Fishing, data=FDenv,ntree=1000,importance=T, mtry=2)
  par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
  partialPlot(fit, x.var=Depth,FDenv,main="")
  partialPlot(fit, x.var=SBT,FDenv,main="")
  partialPlot(fit, x.var=SBS,FDenv,main="")
  partialPlot(fit, x.var=Chl,FDenv,main="")
  partialPlot(fit, x.var=SBT_sea,FDenv,main="")
  partialPlot(fit, x.var=Chl_sea,FDenv,main="")
  partialPlot(fit, x.var=Fishing,FDenv,main="")
  par(mfrow=c(1,1),mar=c(5,4,2,2))
  return(print(fit))
}
funRf(Metrics[9])# Plot random forest response plots and summary stats

# More in detail using trophic level as an example
fitTLrf<-randomForest(CWM.Trophic.level~Depth+SBT+SBS+Chl+SBT_sea+Chl_sea+Fishing, data=FDenv,ntree=1000,importance=T, mtry=2)
plot(fitTLrf) # Check error against number of trees used for training = stable after 200 (hence 100 trees is more than enough) 
print(fitTLrf) # Inspect summary stats

# See partial plots - i.e., model predictions across the range of values of each predictor (while keeping all other covariates fixed at mean levels)
par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
partialPlot(fitTLrf, x.var=Depth,FDenv,main="")
partialPlot(fitTLrf, x.var=SBT,FDenv,main="")
partialPlot(fitTLrf, x.var=SBS,FDenv,main="")
partialPlot(fitTLrf, x.var=Chl,FDenv,main="")
partialPlot(fitTLrf, x.var=SBT_sea,FDenv,main="")
partialPlot(fitTLrf, x.var=Chl_sea,FDenv,main="")
partialPlot(fitTLrf, x.var=Fishing,FDenv,main="")
par(mfrow=c(1,1),mar=c(5,4,2,2))
# Check variable importance plots
varImpPlot(fitTLrf)

#################################################
# END
################################################
# Optional
########################################3
# Test GLMM and extract parameters to fascilitate comparisons across organism groups and areas
# Please not that you may need to introduce square terms to account for non-linearities

library(nlme)
fitTLglm <- lme(CWM.Trophic.level~Depth+SBT+SBS+Chl+SBT_sea+Chl_sea+Fishing,
                random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat), # Can try other error structures
                data = FDenv, control=lmeControl(opt="optim"))
summary(fitTLglm)
anova(fitTLglm)

# Extract parameters
str(summary(fitTLglm))
summary(fitTLglm)$coefficients


##Plot some diagnostics
plot(fitTLglm)
qqnorm(residuals(juvenile))

#observed versus fitted values
plot(fitTLglm, CWM.Trophic.level~ fitted(.), abline = c(0,1))

################################################################
