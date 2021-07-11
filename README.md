# TWAS
Multi-tissue TWAS analysis was performed using transcriptome prediction models from adipose subcutaneous, adipose visceral omentum, ovary, uterus, vagina and whole blood.
All computer code used in this study are published and freely available in the referenced articles and their URLS are:
S-PrediXcan and S-MultiXcan https://github.com/hakyimlab/MetaXcan, 
COLOC https://cran.r-project.org/web/packages/coloc/index.html, 
JTI and MR-JTI https://github.com/gamazonlab/MR-JTI.

# S-PrediXcan analysis
We performed S-PrediXcan analysis using endometrial cancer GWAS summary statistics, tissue transcriptome model and covariance. 
Tissue transcriptome model and covariance was obtained from https://github.com/hakyimlab/MetaXcan. 

```bash
module load python/3.6.1

$ ./SPrediXcan.py \
--model_db_path ${path to tissue transcriptome model} \
--covariance ${path to file containing the covariance matrices for the SNP dosage and transcriptome} \
--gwas_folder data/GWAS \
--gwas_file_pattern "allEC.gz" \ #headers are SNPID,POS,CHR,OA_INFO,MEAN_EAF,EA,OA,BETA,SE,P,N
--snp_column SNPID \
--effect_allele_column EA \
--non_effect_allele_column OA \
--beta_column BETA \
--pvalue_column P \
--output_file results/spredixcan_allEC_${tissue}.csv \
--throw 
```

# S-MultiXcan analysis
We performed S-MultiXcan analysis using S-PrediXcan results from adipose subcutaneous, adipose visceral omentum, ovary, uterus, vagina and whole blood.

```bash
module load python/3.6.1

./SMulTiXcan.py \
--models_folder data/models_v8p/elastic_net_models/ \
--models_name_filter ".*Adipose_Subcutaneous.*db" ".*Adipose_Visceral_Omentum.*db" ".*Ovary.*db" ".*Uterus.*db" ".*Vagina.*db" ".*Whole_Blood.*db" \
--models_name_pattern "en_(.*).db" \
--snp_covariance data/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz \ #gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz was obtained from https://github.com/hakyimlab/MetaXcan
--metaxcan_folder results/ \
--metaxcan_filter "spredixcan_allEC_(.*).csv" \
--metaxcan_file_name_parse_pattern "spredixcan_allEC_(.*).csv" \
--gwas_file data/GWAS/allEC.gz \
--snp_column SNPID --non_effect_allele_column OA --effect_allele_column EA --beta_column BETA --pvalue_column P --se_column SE \
--cutoff_condition_number 30 \
--verbosity 7 \
--throw \
--output results/smultixcan_allEC.csv
```

# Coloc analysis
For MultiXcan-identified genes, we performed COLOC analysis using their cis-eQTL summary statistics from GTEx v8 (https://gtexportal.org/home/) and their association estimates on endometrial cancer GWAS.

```R
module load R/4.0.2

library(coloc)
args <- commandArgs(TRUE)
gwas <- read.table("allEC.tsv", header = T) #headers are snp,a1,a2,maf,beta.gwas,se.gwas,p.gwas
eqtl <- read.table("eQTL.tsv", header = T) #headers are gene,snp,a1,a2,beta.eqtl,se.eqtl,p.eqtl
for(gene in unique(eqtl$gene)){
        print(gene)
        temp <- eqtl[eqtl$gene==gene,]
        gwas_temp <- subset(gwas, snp %in% temp$snp)
        data <- merge(gwas_temp, temp, by = "snp")
	my.res <- coloc.abf(dataset1=list(beta=data$beta.gwas, varbeta=data$se.gwas*data$se.gwas, N=121885,s=0.10588669647,type="cc"), 
		    dataset2=list(beta=data$beta.eqtl, varbeta=data$se.eqtl*data$se.eqtl, N=${sample size}, type = "quant"), 
		    MAF = data$maf)
}
```


# JTI analysis 
We performed JTI analysis using endometrial cancer GWAS, tissue transcriptome model and covariance.
Tissue transcriptome model and covariance were obtained from https://zenodo.org/record/3842289#.YLxPffkzaCo.

```bash

module load python/3.6.1

$ ./SPrediXcan.py \
--model_db_path ${path to tissue transcriptome model} \
--covariance ${path to file containing the covariance matrices for the SNP dosage and transcriptome} \
--gwas_folder data/GWAS \
--gwas_file_pattern "allEC.gz" \
--snp_column SNPID \
--effect_allele_column EA \
--non_effect_allele_column OA \
--beta_column BETA \
--pvalue_column P \
--output_file results/JTI_allEC_${tissue}.csv \
--throw 
```

# MR-JTI analysis
For JTI-identified genes, we performed MR-JTI analysis using their cis-eQTL summary statistics from GTEx v8 (https://gtexportal.org/home/) and their association estimates on endometrial cancer GWAS. (https://gtexportal.org/home/).
All cis-eQTLs for each gene were clumped using Plink (https://www.cog-genomics.org/plink/), before running MR-JTI analysis. 


```bash
module load plink/1.90b6.8

plink \
--bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC \
--clump ${eQTL summary statistics of JTI-identified gene} \ #headers are SNP,CHR,bp,A1,A2,beta,se,P
--clump-r2 0.01 \
--out ${clumped eQTLs} \
--clump-p1 1 \
--thread-num 1

```

```R
module load R/3.4.1

Rscript MR-JTI.r \
--df_path ${dataframe of GWAS and eQTL summary statistics} \ #headers are rsid, effect_allele_gwas, ldscore, eqtl_beta, eqtl_se, eqtl_p, gwas_beta, gwas_se, gwas_p
--n_genes ${total number of genes tested in each tissue} \
--out_path results_${gene}.csv \

```
