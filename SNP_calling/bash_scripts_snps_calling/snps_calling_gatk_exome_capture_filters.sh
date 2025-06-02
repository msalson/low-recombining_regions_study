#!/bin/bash
#SBATCH -J filter
#SBATCH -o filter.%j.out
#SBATCH -e filter.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH -c 10
#SBATCH --mem 10G

module load gatk4/4.2.3.0
module load bcftools/1.10.2
module load vcftools/0.1.16

bcftools concat --file-list all_chr.txt --output All_chr_all_samples.vcf --threads 10

bcftools view --types snps -m 2 -M 2 All_chr_all_samples.vcf -o All_chr_all_samples_onlyBiallelic.vcf

gatk IndexFeatureFile -I All_chr_all_samples_onlyBiallelic.vcf

gatk VariantFiltration --java-options -Xmx10G -R reference.fasta -O All_chr_all_samples_onlyBiallelic_filtered.vcf -V All_chr_all_samples_onlyBiallelic.vcf --cluster-size 3 --cluster-window-size 10 --filter-expression "QUAL<60.0" --filter-name "LOW_QUAL" --filter-expression "DP<650.0" --filter-name "LOW_DP" --filter-expression "DP>13000.0" --filter-name "HIGH-DP" --filter-expression "FS>60.0" --filter-name "HIGH-FS" --filter-expression "QD<2.0" --filter-name "LOW-QD" 

bcftools view -f "PASS" -o All_chr_all_samples_onlyBiallelic_filtered_PASS.vcf All_chr_all_samples_onlyBiallelic_filtered.vcf

vcftools --vcf All_chr_all_samples_onlyBiallelic_filtered_PASS.vcf --minDP 5 --maxDP 100 --recode --recode-INFO-all --out All_chr_all_samples_onlyBiallelic_filtered_PASS_masking5-100

bcftools view -e 'F_MISSING>=0.5' --threads 5 -Ov -o All_chr_all_samples_onlyBiallelic_filtered_PASS_masking5-100_removing05missing.vcf All_chr_all_samples_onlyBiallelic_filtered_PASS_masking5-100.recode.vcf


