# Teaching materials for the GWAS practice  

## Description  
The set of source codes to prepare input files and slides for the GWAS practice at BIO373.  
Here we see how GWAS/GS works with rice agronomic traits.  

To prepare input genotypes and run GWAS/GS, 
1. impute genotypes using BEAGLE (PLINK2VCF2BEAGLE.sh)  
2. covert the imputed vcf into a numeric format of GAPIT input (GAPIT_GenoPrep.R)  
3. run GAPIT using "GAPIT_practice.R"  

To make slides,  
1. clone this repository in your local environment
2. make a RStudio project
3. knit GWAS_practice_Yasu.Rmd as a beamer presentation


## References 
- GAPIT, https://zzlab.net/GAPIT/
- Zhao, Keyan, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, et al. “Genome-Wide Association Mapping Reveals a Rich Genetic Architecture of Complex Traits in Oryza Sativa.” Nature Communications 2, no. 1 (September 2011): 467. https://doi.org/10.1038/ncomms1467.  
- Rice Diversity DB, http://www.ricediversity.org/data/  
