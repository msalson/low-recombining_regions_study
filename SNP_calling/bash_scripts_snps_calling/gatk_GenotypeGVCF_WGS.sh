#!/bin/bash
#SBATCH -J Ggvcf_XX
#SBATCH -o Ggvcf_XX.%j.out
#SBATCH -e Ggvcf_XX.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH -c 1          
#SBATCH --mem 60G

module load gatk4/4.2.3.0

gatk --java-options -Xmx55G GenotypeGVCFs -R pearl_millet_ONT.fasta -V gendb://mydbSuper-Scaffold_XX -O Super-Scaffold_XX_173samples.vcf.gz

