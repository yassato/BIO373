# This is a shell script to prepare imputed SNP data
# Install PLINK and BEAGLE first
# For instalation, bioconda channel of the anaconda can used as
# "conda install -c bioconda plink or beagle"

# Data are downloaded from　http://www.ricediversity.org/data/sets/44kgwas/
# Original data from Zhao et al. (2011) Nat. Comm.
cd RiceDiversity_44K_Genotypes_PLINK

# Convert plink to vcf
plink --file sativas413 --recode vcf

# SNP imputation using BEAGLE
beagle gt=plink.vcf out=RiceDiversity_44K_Genotypes_PLINK_imputed
mv RiceDiversity_44K_Genotypes_PLINK_imputed.vcf.gz ../

cd ../
mv RiceDiversity_44K_Genotypes_PLINK_imputed.vcf.gz ./
