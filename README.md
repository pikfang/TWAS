# testing
```bash
##MetaXcan analysis##

module load python/3.6.1

#Tissue transcriptome model and covariance was obtained from https://github.com/hakyimlab/MetaXcan 
#We performed PrediXcan analysis using transcriptome prediction models from adipose subcutaneous, adipose visceral omentum, ovary, uterus, vagina and whole blood.
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
--output_file results/spredixcan_allEC_${tissue}.csv
--throw 

##MultiXcan analysis##

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

```R
##Coloc analysis##
module load R/4.0.2

library(coloc)
args <- commandArgs(TRUE)
gwas <- read.table(args[1], header = T) #headers are snp,a1,a2,beta.gwas,se.gwas,p.gwas
eqtl <- read.table(args[2], header = T) #headers are gene,snp,a1,a2,beta.eqtl,se.eqtl,p.eqtl
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

```bash
##JTI analysis##

#Tissue transcriptome model and covariance was obtained from https://zenodo.org/record/3842289#.YLxPffkzaCo
#We performed JTI analysis using transcriptome prediction models from adipose subcutaneous, adipose visceral omentum, ovary, uterus, vagina and whole blood.

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
--output_file results/JTI_allEC_${tissue}.csv
--throw 
```
