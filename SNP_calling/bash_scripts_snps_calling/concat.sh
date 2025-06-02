#!/bin/bash
#SBATCH -J concat
#SBATCH -o concat.%j.out
#SBATCH -e concat.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH -c 10

module load bcftools/1.15.1

bcftools concat --file-list names_files_scaffolds.txt --output ALL_scaffolds_173_samples.vcf --threads 10
