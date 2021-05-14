##############################
#GWAS/GS practice using GAPIT#
##############################

# download GAPIT source code
source("https://zzlab.net/GAPIT/gapit_functions.txt")

# download phenotype data
p <- read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
                sep="\t", header=TRUE)

# read genotype data
g <- read.table("./data/RiceDiversity_44K_Genotypes_PLINK_imputed.txt.gz",
                header=TRUE,sep="\t")
gm <- read.table("./data/RiceDiversity_44K_Genotypes_PLINK_info.txt"
                 ,header=TRUE,sep="\t")

# GWAS with a naive linear model
myGAPIT_GLM <- GAPIT(
  Y=p[,c("HybID","Seed.length")],
  GD=g,
  GM=gm,
  Inter.Plot=TRUE,
  model="GLM")

# GWAS with a linear mixed model
myGAPIT_MLM <- GAPIT(
  Y=p[,c("HybID","Seed.length")],
  GD=g,
  GM=gm,
  Inter.Plot=TRUE,
  model="MLM")


# gBLUP for the flowering time 2006 at Arkansas
myGAPIT_BLUP <- GAPIT(
  Y=p[,c("HybID","Year06Flowering.time.at.Arkansas")],
  GD=g,
  GM=gm,
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
plot(pred$Prediction,
     y$Year07Flowering.time.at.Arkansas,
     ylab="flowering on 2007", xlab="predicted")

cor.test(pred$Prediction,
         y$Year07Flowering.time.at.Arkansas)

# Predicting flowering time of missing accessions
plot(pred$Prediction[is.na(y$Year06Flowering.time.at.Arkansas)],
     y$Year07Flowering.time.at.Arkansas[is.na(y$Year06Flowering.time.at.Arkansas)],
     ylab="flowering on 2007", xlab="predicted")

cor.test(pred$Prediction[is.na(y$Year06Flowering.time.at.Arkansas)],
         y$Year07Flowering.time.at.Arkansas[is.na(y$Year06Flowering.time.at.Arkansas)])


