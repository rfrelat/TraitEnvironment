---
title: "Introduction to RLQ analysis"
author: "Frelat, R."
date: "7th June 2021"
output:
  html_document: 
    keep_md: yes
  word_document: default
  pdf_document: default
---



This document provides an introduction to RLQ analysis. The tutorial targets students and scientists in marine biology and ecology with previous knowledge of the [R software](https://cran.r-project.org/). 

Please also consult the tutorial by Dray S. XXX for a detailed on RLQ analysis. 

# 1. Preliminaries

## 1.1. Load ade4 package

The analyses require the R packages [ade4 (v ≥ 1.7.16)](https://pbil.univ-lyon1.fr/ade4/home.php?lang=eng).


```r
library(ade4)
```

If you want to plot a map with country border, then you would also need the packages `maps (v ≥ 3.3)` and `mapdata (v ≥ 2.3)`.


```r
library(maps); library(mapdata)
```

If you get an error message, check that the R packages are installed correctly. If not, use the command: `install.packages(c("ade4", "maps", "mapdata"))`.

The example dataset is available as the Rdata file `NorthSea_FishTraitEnv.Rdata`, available for download [here](https://github.com/rfrelat/TraitEnvironment/raw/main/NorthSea_FishTraitEnv.Rdata).  
(https://github.com/rfrelat/TraitEnvironment/raw/main/RLQanalysis.html).  

## 1.2 Load the example dataset

Make sure the file `NorthSea_FishTraitEnv.Rdata` is in your working directory, then load it in R.


```r
load("NorthSea_FishTraitEnv.Rdata")
```

The Rdata file contains four objects: 
- `abu` containing the abundance of taxa in grid cells
- `env` containing the environmental condition per grid cell
- `trait` containing the trait information per taxa 
- `coo`: the coordinates of each grid cell


Importantly, the rows in `abu` correspond to the same grid cell than the rows in `env`, and the column in `abu` correspond to the same taxa than the rows in `trait`.  
If you want to learn how to create such dataset, see the short tutorial on setting trait-environement dataset in XXX.


```r
all(row.names(abu)==row.names(env))
```

```
## [1] TRUE
```

```r
all(colnames(abu)==row.names(trait))
```

```
## [1] TRUE
```

Using this fish community of the North Sea as an example, you will learn how to compute the RLQ analysis.

## 1.3 Quick summary of the variables


```r
dim(trait)
```

```
## [1] 90  7
```

```r
names(trait)
```

```
## [1] "Trophic level"  "K"              "Lmax"           "Lifespan"      
## [5] "Offspring size" "Fecundity"      "Age maturity"
```

Seven traits for 90 taxa which are ...



```r
dim(env)
```

```
## [1] 111   7
```

```r
names(env)
```

```
## [1] "Depth"   "SBT"     "SBS"     "Chl"     "SBT_sea" "Chl_sea" "Fishing"
```

Seven environemental variables for 111 grid cells which are ...

# 2. Single multivariate analysis

## 2.1 COA on abundance matrix L (sites x species)

Some explanation on COA and what is the difference between PCA and COA


```r
dim(abu) # 111 sites x 90 species
```

```
## [1] 111  90
```

```r
coa.abu <- dudi.coa(abu, scannf = FALSE, nf=2)
```

## 2.2 PCA on traits matrix Q (species x traits)

Using species weight from previous COA 

Mention that it could be mix dataset (Hill and Smith analysis)


```r
dim(trait) # 90 species x 7 continuous traits
```

```
## [1] 90  7
```

```r
pca.trait<-dudi.pca(trait, scannf = FALSE, 
                     row.w = coa.abu$cw) 
```

## 2.3 PCA on environment matrix R (sites x environment)

using sites weight from previous COA 

Same here, it could be categorical variables
And the number of environmental variables do not have to match the number of traits (indeed, easier if it doesn't match)


```r
dim(env) # 111 sites x 7  environmental variables
```

```
## [1] 111   7
```

```r
pca.env <- dudi.pca(env, scannf = FALSE, 
                     row.w = coa.abu$lw)
```


# 3. Run RLQ analysis and interpret the results

## 3.1 Compute RLQ analysis


```r
rlqF <- rlq(pca.env, coa.abu, pca.trait, 
               scannf = FALSE)
```



## 3.2 Visualize the output


```r
#Plot traits score
t1 <- order(rlqF$c1[,1])
dotchart(rlqF$c1[t1,1], pch=16, 
         labels = names(trait)[t1])
abline(v=0, lty=2)
```

![](RLQanalysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


```r
#Plot environment score
e1 <- order(rlqF$l1[,1])
dotchart(rlqF$l1[e1,1], pch=16,
         labels = names(env)[e1])
abline(v=0, lty=2)
```

![](RLQanalysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


## 3.2 Visualize the output


```r
#Choice of diverging color scale (from RColorBrewer)
colpal <- terrain.colors(6)
#Choice of diverging color scale (from RColorBrewer)
# colpal <- rev(brewer.pal(6,"RdYlBu")) 
#Set the color of each site based on RLQ score
colo <- colscale(rlqF$lR[,1], colpal)

#Plot the map
apply(coo, 2, range)
```

```
##      Longitude Latitude
## [1,]      -5.5     51.5
## [2,]      12.5     61.5
```

```r
map("world", xlim=c(-6, 13), ylim=c(50, 63), 
    mar=c(5,5,1,0), col="gray70", fill=TRUE, 
    border=NA, lforce="e")
points(coo, col=colo$col, pch=".", cex=6)
map.axes()
```

![](RLQanalysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


#4. Fourth corner analysis

Fourth corner analysis takes very long to compute.
Additionally, multiple repetitions are needed in order to estimate the p-values with enough precision before p-values are adjusted for multiple testing.

*might not be usefull here...*


```r
nrep <- 1000
fc <- fourthcorner(env, abu, trait, nrepet=nrep)
plot(fc)
```

![](RLQanalysis_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


# References

Dray, S., Choler, P., Dolédec, S., Peres-Neto, P. R., Thuiller, W., Pavoine, S., & ter Braak, C. J. (2014). Combining the fourth‐corner and the RLQ methods for assessing trait responses to environmental variation. Ecology, 95(1), 14-21. [DOI 10.1890/13-0196.1](https://doi.org/10.1890/13-0196.1) 