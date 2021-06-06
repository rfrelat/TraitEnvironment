############################################################
# FutureMares Task 1.2. Trait-environment workshop. 7th June 2021
# Example data and R code for estimating HSMC 
# Benjamin Weigel, 7/6/2021.
##############################################################

knitr::opts_chunk$set(echo = TRUE)

library(Hmsc) 
library(tidyverse)
library(viridis)
library(vioplot)
library(abind)
library(RColorBrewer)
library(ape)
library(corrplot)

load("NEAtl_FishTraitEnv.Rdata")

all(row.names(abu)==row.names(env))
all(colnames(abu)==row.names(trait))

dim(trait)
names(trait)

trait<-trait %>%  
  dplyr::rename_all(make.names)
head(trait)

dim(env)
names(env)

str(taxo) #listed as chr
taxo <- as.data.frame(unclass(taxo)) # make factors
str(taxo) # now they are!

# build the tree  
tree<-as.phylo(~class/order/family/genus/species, data=taxo, collapse = FALSE) 
tree$edge.length<-rep(1,length(tree$edge))

# It's important to check if the tip lables of tree correspond to the names in abu
if(all(sort(tree$tip.label) == sort(colnames(abu)))){
  print("species names in tree and abu match")
} else{
  print("species names in tree and abu do not match")
} #there is some problem which we can already see from the NA message before.

# Looks like we have some NA in genus, so what we do now is take the species information for all NA at genus level to fill the NA correctly.
taxo2<-taxo%>% 
  mutate(genus = coalesce(genus,species))

# Build tree again with new taxo2
tree<-as.phylo(~class/order/family/genus/species, data=taxo2, collapse = FALSE) # no NA message
tree$edge.length<-rep(1,length(tree$edge))

# Check lable correspondence
 if(all(sort(tree$tip.label) == sort(colnames(abu)))){
  print("species names in tree and abu match")
} else{
  print("species names in tree and abu do not match")
} # All good now!

#Have a look at the tree
str(tree)
plot(tree, cex=0.5)

XFormula = as.formula(paste("~",paste(colnames(env), collapse="+")))
XFormula

TrFormula = as.formula(paste("~",paste(colnames(trait), collapse="+")))
TrFormula

studyDesign = data.frame(grid.cell = as.factor(rownames(coo))) # "sample" level is defined as grid.cell
rL = HmscRandomLevel(sData = coo) # here the random effect is defined to be spatially explicit based on the coordinates of the grid cell
rL = setPriors(rL,nfMin=1,nfMax=2) # here we set the number of latent variable factors. In this example they are low for computational reasons 

# Have a look at the random effect 
rL
head(rL$s)


test.run = TRUE # TEST RUN TAKES ABOUT 2 MIN ON MY COMPUTER
if (test.run){
  # with this option, the model runs fast but results are not reliable
  thin = 1 # interval of iterations per samples 
  samples = 10 # number of samples taken per chain
  transient = ceiling(0.5*samples*thin) # burn-in, i.e. number of first iterations which will be cut off
  nChains = 2 # number of independent MCMC chains 
} else { 
  # with this option, the model evaluates slow but it reproduces the results shown in this tutorial
  thin = 100 # interval of iterations per samples 
  samples = 250 # number of samples taken per chain
  transient = ceiling(0.5*samples*thin) # burn-in, i.e. number of first iterations which will be cut off
  nChains = 4 # number of independent MCMC chains 
  }

m = Hmsc(Y= abu, XData = env,  XFormula = XFormula, TrFormula = TrFormula, TrData = trait, phyloTree = tree, studyDesign = studyDesign, ranLevels = list("grid.cell"= rL),  distr = "normal")

##  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains)
## 
## # You can also run chains in parallel by adding 'nParallel = nChains' in the sampleMCMC() function, but you won't see the progress of the run until it finished.
## 
## filename =  paste("NorthAtlantic","_thin_", as.character(thin),"_samples_", as.character(samples),"_chains_", as.character(nChains), ".Rdata",sep = "")
## 
## save(m,file=filename)

load("NorthAtlantic_rL_tree_thin_10_samples_250_chains_4.Rdata")

# convert model to coda object
mpost = convertToCodaObject(m)

# check effective sample size 
# should ideally be close to your drawn samples, here, 4 (chains) * 250 (samples) = 1000 samples)
summary(effectiveSize(mpost$Beta)) 

# check the distribution of the effective sample size (ess) and the Gelman-Rubin potential scale reduction factor (psrf), the latter should be ideally (very) close to 1, if not, chains did not properly converge and need more iterations. 
# These should be checked for all parameters of interest, i.e. beta and gamma
# Provided model looks quite okay already but could be better.

par(mfrow=c(1,3))
# Beta parameters (species-environment) 
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf, main="psrf(beta)")
vioplot(gelman.diag(mpost$Beta,multivariate = FALSE)$psrf,main="psrf(beta)")

# Gamma parameters (trait-environment)
hist(effectiveSize(mpost$Gamma), main="ess(gamma)")
hist(gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf, main="psrf(gamma)")
vioplot(gelman.diag(mpost$Gamma,multivariate = FALSE)$psrf,main="psrf(gamma)")

# compute predicted species abundance matrix/ posterior samples 
predY=computePredictedValues(m)

MF = evaluateModelFit(hM = m, predY = predY)
MF
mean(MF$R2)

#Here, we take the mean of posterior samples for further processing
predY = apply(abind(predY,along=3),c(1,2), mean)

# Plot species specific R2 in relation to prevalence
plot(colSums(((m$Y>0)*1)/m$ny), MF$R2,main=paste("Mean R2 = ", round(mean(MF$R2),2),".", sep=""), xlab = "Prevalence")

par(mfrow=c(1,1))
# Specify groups of how the variation should be partitioned You can also combine groups 
group=c(1,2,3,4,5,6,7)
# Specify group names m$covNames[-1] gives the included covariate names excluding the intercept
groupnames = m$covNames[-1]
#compte species specific variance partitioning
VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)

plotVariancePartitioning(m, VP, viridis(8)) 

# Shown are all species with 90% posterior support for having positive (red) or negative(blue) responses to environmental covariates
beta = getPostEstimate(m, "Beta")
plotBeta(m, beta, supportLevel=.9, spNamesNumbers = c(FALSE, FALSE), covNamesNumbers = c(TRUE, FALSE), plotTree = T)

gamma = getPostEstimate(m, "Gamma")
plotGamma(m, gamma, supportLevel=.8) 
# modify support level to highlight stronger or weaker posterior support for positive/ negative relationship

VP$R2T$Y
# The traits explain only 1.93% of variation in species abundances 

VP$R2T$Beta
barplot(VP$R2T$Beta[-1])
#The traits explain also not much of the variation in the species niches, with traits explaining most out of the variation in species responses to SBT, SBT_sea and fishing with ca. 2.5% each

# Construct environmental Gradient based on fitted model, specify your focal variable. 
Gradient = constructGradient(m, focalVariable="SBT", ngrid = 25)
# make predictions based on fitted model 
predYgradient = predict(m,XData=Gradient$XDataNew, ranLevels = Gradient$rLNew, studyDesign = Gradient$studyDesignNew, expected = TRUE) 

par(mfrow=c(2,4))
# measure = T implies trait relationship to focal variable, index = 2 implies second trait of trait matrix (1 is the intercept). 

# traits ~ SBT 
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 2, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 3, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 4, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 5, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 6, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 7, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 8, showData =F)
# Species  ~ SBT // If you specify measure = Y, this is the species response, with index being the number of species in the abu matrix
plotGradient(m, Gradient, pred=predYgradient, measure="Y", index = 50, showData =F) #here examplified for species 50, cod

# Note that the scale of all relationships is very small for the traits

S=rowSums(predY)
#predicted CWM
predT = (predY%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)

#make data frame for plotting
data<-data.frame(predT, coo)
colpal<-rev(brewer.pal(11,"RdYlBu")) # Color palette for the map

funPlot<-function(Var){
colnames(data)[which(colnames(data)==Var)]<-"Var"
ggplot(data, aes(x= Longitude, y= Latitude, fill= Var))+
  geom_tile(data=data,aes(x=Longitude,y=Latitude,fill= Var)) + # 
  scale_fill_gradientn(name = Var, colours=colpal, na.value = 'white')+
  borders(fill="gray44",colour="black") +
  coord_quickmap(xlim=c(range(data$Longitude)),ylim=c(range(data$Latitude)))+
  labs(x = "Lon",y="Lat")+
  theme(legend.position = c(0.93,0.8))+
  theme(legend.title = element_text(size = 8),legend.text = element_text(size = 7))+
  guides(shape = guide_legend(override.aes = list(size = 0.2)))+
  theme(panel.background = element_rect(fill=alpha('light blue', 0.4), colour = 'black'))
}
# list all output 
Metrics<-colnames(data)
Metrics

#We plot the CWM of Trophic level, which is the 2nd variable
funPlot(Var=Metrics[2])# Plot CWM traits and diversity metrics

par(mfrow=c(1,1))
OmegaCor = computeAssociations(m)
supportLevel = 0.9
for (r in 1:m$nr){ 
  plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE") 
  toPlot = ((OmegaCor[[r]]$support>supportLevel) +                                    (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean 
  par(xpd=T) 
  colnames(toPlot)=rownames(toPlot)=gsub("_"," ",x=colnames(toPlot)) 
  corrplot(toPlot[plotOrder,plotOrder], method = "color", col=colorRampPalette(c("blue","white","red"))(200), title=paste("random effect level:",m$rLNames[r]),type="full",tl.col="black",tl.cex=.4, mar=c(0,0,6,0))
}  

m = Hmsc(Y= 1*(abu>0), XData = env,  XFormula = XFormula, TrFormula = TrFormula, TrData = trait, phyloTree = tree, studyDesign = studyDesign, ranLevels = list("grid.cell"= rL),  distr = "probit")

##  m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains)

load("NorthAtlantic_rL_tree_PA_thin_10_samples_250_chains_4.Rdata")

# compute predicted species occurrence matrix from posterior samples 
predY=computePredictedValues(m)

#evaluate model fit
MF = evaluateModelFit(hM = m, predY = predY)
MF
mean(MF$AUC)
mean(MF$TjurR2)

#Here, we take the mean of posterior samples for further processing
predY = apply(abind(predY,along=3),c(1,2), mean)

# Plot species specific R2 in relation to prevalence
plot(colSums(((m$Y>0)*1)/m$ny), MF$AUC,main=paste("Mean AUC = ", round(mean(MF$AUC),2),".", sep=""), xlab = "Prevalence")

par(mfrow=c(1,1))
# Specify groups of how the variation should be partitioned You can also combine groups 
group=c(1,2,3,4,5,6,7)
# Specify group names m$covNames[-1] gives the included covariate names excluding the intercept
groupnames = m$covNames[-1]
#compte species specific variance partitioning
VP = computeVariancePartitioning(hM = m, group = group, groupnames = groupnames)

plotVariancePartitioning(m, VP, viridis(8)) 

# Shown are all species with 90% posterior support for having positive (red) or negative(blue) responses to environmental covariates
beta = getPostEstimate(m, "Beta")
plotBeta(m, beta, supportLevel=.9, spNamesNumbers = c(FALSE, FALSE), covNamesNumbers = c(TRUE, FALSE), plotTree = T)

gamma = getPostEstimate(m, "Gamma")
plotGamma(m, gamma, supportLevel=.9) 

VP$R2T$Beta

VP$R2T$Y

# Construct environmental Gradient based on fitted model, specify your focal variable. 
Gradient = constructGradient(m, focalVariable="SBT", ngrid = 25)
# make predictions based on fitted model 
predYgradient = predict(m,XData=Gradient$XDataNew, ranLevels = Gradient$rLNew, studyDesign = Gradient$studyDesignNew, expected = TRUE) 

par(mfrow=c(2,4))
# measure = T implies trait relationship to focal variable, index = 2 implies second trait of trait matrix (1 is the intercept). 

# traits ~ SBT 
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 2, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 3, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 4, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 5, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 6, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 7, showData =F)
plotGradient(m, Gradient, pred=predYgradient, measure="T", index = 8, showData =F)

#in a presence absence probit model, measure S provides the response of species richness against your chosen covariate
plotGradient(m, Gradient, pred=predYgradient, measure="S", index = 1, showData =F) # 
