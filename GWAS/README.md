# Teaching materials for the GWAS practice  

## Description  
The set of source codes to prepare input files and slides for BIO373 GWAS practice.  
To prepare input genotypes and run GWAS/GS of rice agronomic traits, 
1. Impute genotypes using BEAGLE (./preprocessing/PLINK2VCF2BEAGLE)  
2. Covert the imputed vcf into a numeric format of GAPIT input (.data/GAPIT_GenoPrep.R)  
3. Run GAPIT using "GAPIT_practice.R"    


## References 
- GAPIT, https://zzlab.net/GAPIT/
- Zhao, Keyan, Chih-Wei Tung, Georgia C. Eizenga, Mark H. Wright, M. Liakat Ali, Adam H. Price, Gareth J. Norton, et al. “Genome-Wide Association Mapping Reveals a Rich Genetic Architecture of Complex Traits in Oryza Sativa.” Nature Communications 2, no. 1 (September 2011): 467. https://doi.org/10.1038/ncomms1467.  
- Rice Diversity DB, http://www.ricediversity.org/data/  
