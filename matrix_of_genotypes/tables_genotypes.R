#!/usr/bin/env Rscript

# Script to perform the analysis to plot the matrix of genotypes.
# For the PCAdapt analysis, the vcf file had to be converted in a bed format with plink (v. 1.9) and --make-bed parameter: plink --vcf spatial_dataset_chr3_101-189Mb.vcf --make-bed --out spatial_dataset_chr3_101-189Mb_BED

library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(scales) #Version: 1.2.1
library(pcadapt) #Version: 4.3.3

#1. PCAdapt analysis to obtain the SNPs contributed the most to the population structure:
bed <- "data/spatial_dataset_chr3_101-189Mb_BED.bed"
acp <- read.pcadapt(bed, type='bed')
pcadapt <- pcadapt(acp, K = 2)

# get the outliers
padj <- p.adjust(pcadapt$pvalues, method = "BH")
alpha <- 0.05
outliersBH0.05 <- which(padj < alpha)
snp_pc <- get.pc(pcadapt, outliersBH0.05)
colnames(snp_pc) <-  c("SNP_number", "PC")
snp_pc <- snp_pc %>% rowwise() %>%
  mutate(sequence = colnames(acp)[SNP_number])

write.table(snp_pc,file="SNPs_outliers_PCADAPT.txt",col.names=FALSE)
outlierSNPs<- "SNPs_outliers_PCADAPT.txt"

# 2. A VCF file with the samples ordered following the clusters identified in the population has to be used: 
vcfOrdered<-"data/spatial_dataset_chr3_101-189Mb_ORDERED.vcf"

# 3. launching the python script to generate the table of genotypes:
system(paste0("./generate_tables_of_genotypes.py -v", vcfOrdered, " -s ", outlierSNPs, " -m matrix -p matrix -f matrix"))
system("rm matrix")
