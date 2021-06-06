############################################################
# FutureMares Task 1.2. Trait-environment workshop. 7th June 2021
# Example data and R code for RLQ analysis
# Romain Frelat and Esther Beukhof, 7/6/2021.
##############################################################

library(ade4)

library(ggplot2)

# Load data
load("NEAtl_FishTraitEnv.Rdata")

# Check is grid cells match and species names match
all(row.names(abu)==row.names(env))
all(colnames(abu)==row.names(trait))

# Check trait data
dim(trait)
names(trait)

# Check environmental data
dim(env)
names(env)

# 1. Correspondence analysis on L
coa.abu <- dudi.coa(abu, scannf = FALSE, nf=2)

# 2. Multivariate analysis (here, PCA) on Q
pca.trait <- dudi.pca(trait, scannf = FALSE, 
                      row.w = coa.abu$cw)

# 3. Multivariate analysis (here, PCA) on R
pca.env <- dudi.pca(env, scannf = FALSE, 
                     row.w = coa.abu$lw)

# 4. Run RLQ analysis, which compares the co-variance of the three previous steps using co-inertia analysis
rlqF <- rlq(pca.env, coa.abu, pca.trait, 
            scannf = FALSE)

summary(rlqF)

#Plot traits score
t1 <- order(rlqF$c1[,1])
dotchart(rlqF$c1[t1,1], pch=16, 
         labels = names(trait)[t1])
abline(v=0, lty=2)

# top 10 species with positive score
rlqF$mQ[order(rlqF$mQ[,1], decreasing = TRUE)[1:10],]

# top 10 species with negative score
rlqF$mQ[order(rlqF$mQ[,1])[1:10],]

#Plot environment score
e1 <- order(rlqF$l1[,1])
dotchart(rlqF$l1[e1,1], pch=16,
         labels = names(env)[e1])
abline(v=0, lty=2)

# Choice of diverging color scale
colpal <- terrain.colors(7)[-7]
# or from RColorBrewer package
# colpal <- rev(RColorBrewer::brewer.pal(6,"RdYlBu")) 

# Plot the sites scores
mapggplot(coo[,1], coo[,2], rlqF$lR[,1], 
          colpal, main="RLQ1")
