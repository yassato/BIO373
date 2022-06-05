# GWAS/GS exercise using GAPIT
Yasuhiro Sato on dd mm 2022

## 1. Setup and load packages

Clean up a workplace first.
```
rm(list=ls())
```

Load GAPIT source code and its dependency. Wait ca. 5 min. to install everything.  
Some packages are not installed but they are negligible.  
```
source("gapit_functions.txt")
```

## 2. Load phenotypes and genotypes

Download phenotype data
```
pheno_url <- "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"
p <- read.table(pheno_url, sep="\t", header=TRUE)
nrow(p) # No. of plants
head(p)
```

Read genotype data and marker information
```
g <- read.table("RiceDiversity_44K_Genotypes_PLINK_imputed.txt.gz",
                header=TRUE, sep="\t")
gm <- read.table("RiceDiversity_44K_Genotypes_PLINK_info.txt.gz",
                 header=TRUE, sep="\t")
nrow(g[,-1]) # No. of SNPs
ncol(g) # No. of plants
head(gm) # marker info
```

## 3. Run GWAS

Compare the general linear model (GLM) and mixed linear model (MLM).  
Some warmings occur but the analysis still works. When they are finished, output files appear in the current directory.  
```
myGAPIT <- GAPIT(
  Y=p[,c("HybID","Seed.length")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05, # cut-off minor alleles at 0.05
  Inter.Plot=TRUE,  # option to make interactive plots
  model=c("GLM", "MLM"),
  Multiple_analysis=TRUE)
```

## 4. gBLUP

Calculate BLUP for the flowering time 2006 at Arkansas.  
Some warnings occur but the analysis still works.  
```
myGAPIT_BLUP <- GAPIT(
  Y=p[,c("HybID","Year06Flowering.time.at.Arkansas")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05,
  model="gBLUP")
```

Load results of genomic prediction  
```
pred <- read.csv("GAPIT.MLM.Pred.result.csv")
head(pred)
```

Align predicted and observed traits following the taxa name  
```
pred <- pred[order(pred$Taxa),]
y <- p[order(p$HybID),]
```

Pearson's correlation between predicted and observed flowering  
```
cor.test(pred$Prediction, y$Year06Flowering.time.at.Arkansas)
```

Predicting flowering time 2007  
```
res <- lm(y$Year07Flowering.time.at.Arkansas~pred$Prediction)
plot(pred$Prediction,
     y$Year07Flowering.time.at.Arkansas,
     ylab="flowering 2007", xlab="predicted",
     main=paste("r =",round(sqrt(summary(res)$r.squared),2)))
abline(res)
```

Predicting flowering time of missing accessions
```
NA06 <- is.na(y$Year06Flowering.time.at.Arkansas)
res <- lm(y$Year07Flowering.time.at.Arkansas[NA06]~pred$Prediction[NA06])
plot(pred$Prediction[NA06],
     y$Year07Flowering.time.at.Arkansas[NA06],
     ylab="flowering 2007", xlab="predicted",
     main=paste("r =",round(sqrt(summary(res)$r.squared),2)))
abline(res)
```
