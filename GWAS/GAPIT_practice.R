##############################
#GWAS/GS practice using GAPIT#
##############################

# clean up your workplace
rm(list=ls())

# download GAPIT source code
# if errors happen, try twice to load
source("https://zzlab.net/GAPIT/gapit_functions.txt")

# download phenotype data
pheno_url <- "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"
p <- read.table(pheno_url, sep="\t", header=TRUE)
nrow(p) # No. of plants
head(p)

# read genotype data and marker information
g <- read.table("RiceDiversity_44K_Genotypes_PLINK_imputed.txt.gz",
                header=TRUE, sep="\t")
gm <- read.table("RiceDiversity_44K_Genotypes_PLINK_info.txt.gz",
                 header=TRUE, sep="\t")
nrow(g[,-1]) # No. of SNPs
ncol(g) # No. of plants
head(gm) # marker info

# Run GWAS with a general linear model (GLM) or mixed linear model (MLM)
myGAPIT <- GAPIT( # warnings occur but it still works
  Y=p[,c("HybID","Seed.length")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05, # cut-off minor alleles at 0.05
  Inter.Plot=TRUE,  # option to make interactive plots
  model=c("GLM", "MLM"),
  Multiple_analysis=TRUE)

# When finished, output files appear in the current directory


# gBLUP for the flowering time 2006 at Arkansas
myGAPIT_BLUP <- GAPIT( # warnings occur but it still works
  Y=p[,c("HybID","Year06Flowering.time.at.Arkansas")],
  GD=g,
  GM=gm,
  SNP.MAF=0.05,
  model="gBLUP")

# load results of genomic prediction
pred <- read.csv("GAPIT.MLM.Pred.result.csv")
head(pred)

# align predicted and observed traits following the taxa name
pred <- pred[order(pred$Taxa),]
y <- p[order(p$HybID),]

# Pearson's correlation between predicted and observed flowering
cor.test(pred$Prediction, y$Year06Flowering.time.at.Arkansas)

# Predicting flowering time 2007
res <- lm(y$Year07Flowering.time.at.Arkansas~pred$Prediction)
plot(pred$Prediction,
     y$Year07Flowering.time.at.Arkansas,
     ylab="flowering 2007", xlab="predicted",
     main=paste("r =",round(sqrt(summary(res)$r.squared),2)))
abline(res)

# Predicting flowering time of missing accessions
NA06 <- is.na(y$Year06Flowering.time.at.Arkansas)
res <- lm(y$Year07Flowering.time.at.Arkansas[NA06]~pred$Prediction[NA06])
plot(pred$Prediction[NA06],
     y$Year07Flowering.time.at.Arkansas[NA06],
     ylab="flowering 2007", xlab="predicted",
     main=paste("r =",round(sqrt(summary(res)$r.squared),2)))
abline(res)
