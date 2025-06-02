#!/bin/bash
#SBATCH -J Ggvcf
#SBATCH -o Ggvcf.%j.out
#SBATCH -e Ggvcf.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH -c 1          
#SBATCH --mem 80G

module load gatk4/4.2.3.0

gatk --java-options -Xmx80G GenotypeGVCFs -R ./reference.fasta -V gendb://db_chr1 -O chr1.vcf.gz --include-non-variant-sites

#--include-non-variant-sites

