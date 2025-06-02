#!/bin/bash
#SBATCH -J filt
#SBATCH -o filt.%j.out
#SBATCH -e filt.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem 100G

module load gatk4/4.2.3.0
module load bcftools/1.10.2
module load vcftools/0.1.16

bcftools view --types snps -m 2 -M 2 ALL_scaffolds_173_samples.vcf -o ALL_scaffolds_173_samples_onlyBiallelic.vcf

gatk IndexFeatureFile -I ALL_scaffolds_173_samples_onlyBiallelic.vcf

gatk VariantFiltration --java-options -Xmx100G -R pearl_millet_ONT.fasta -O ALL_scaffolds_173_samples_onlyBiallelic_filtered.vcf -V ALL_scaffolds_173_samples_onlyBiallelic.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression "QUAL<60.0" --filter-name "LOW_QUAL" --filter-expression "DP<865.0" --filter-name "LOW_DP" --filter-expression "DP>17300.0" --filter-name "HIGH-DP" --filter-expression "FS>60.0" --filter-name "HIGH-FS" --filter-expression "QD<2.0" --filter-name "LOW-QD" 

bcftools view -f "PASS" -o ALL_scaffolds_173_samples_onlyBiallelic_filtered_PASS.vcf ALL_scaffolds_173_samples_onlyBiallelic_filtered.vcf

vcftools --vcf ALL_scaffolds_173_samples_onlyBiallelic_filtered_PASS.vcf --minDP 5 --maxDP 100 --recode --recode-INFO-all --out ALL_scaffolds_173_samples_onlyBiallelic_filtered_PASS_masking5-100

bcftools view -e 'F_MISSING>=0.1' --threads 5 -Ov -o ALL_scaffolds_173_samples_onlyBiallelic_filtered_PASS_masking5-100_removing01missing.vcf ALL_scaffolds_173_samples_onlyBiallelic_filtered_PASS_masking5-100.recode.vcf






