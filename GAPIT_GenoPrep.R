####################################
# prep. for genotype data for GAPIT#
####################################

# load imputed SNPs using vcfR package
library(vcfR)
g <- read.vcfR("RiceDiversity_44K_Genotypes_PLINK_imputed.vcf.gz")

# genotype marker information
gm <- as.data.frame(g@fix[,c(3,1,2)])
gm$CHROM <- as.numeric(gm$CHROM)
gm$POS <- as.numeric(gm$POS)

# prepare SNP matrix as a numeric format
g <- g@gt
g[g=="0|0"] = 0 # LM as 0
g[g=="1|0"] = 1 # hetero as 1
g[g=="0|1"] = 1 # hetero as 1
g[g=="1|1"] = 2 # another HM as 2
g <- g[,-1]
g <- as.data.frame(g)
g <- apply(g,2,as.numeric)

# add taxa info
pheno_url <- "http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt"
p <- read.table(pheno_url, sep="\t", header=TRUE)
taxa <- p$HybID
snp_data <- data.frame(taxa,t(g))
rownames(snp_data) <- snp_data$taxa
colnames(snp_data) <- c("taxa",colnames(g))

# write SNP matrix and marker info
write.table(snp_data,"RiceDiversity_44K_Genotypes_PLINK_imputed.txt",sep="\t",row.names=FALSE)
write.table(gm,"RiceDiversity_44K_Genotypes_PLINK_info.txt",sep="\t",row.names=FALSE)

