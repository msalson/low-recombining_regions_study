#!/bin/bash
#SBATCH -J snp
#SBATCH -o snp.%j.out
#SBATCH -e snp.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH -c 1
#SBATCH --mem 100G

module load python/3.9

./snps_from_scaffolds_chr.py -v ALL_scaffolds_173_samples.vcf -w scaffolds_within_chromosomes.csv -o ALL_CHR_173_samples.vcf

sort -k1,2 -V ALL_CHR_173_samples.vcf > ALL_CHR_173_samples_SORTED.vcf
