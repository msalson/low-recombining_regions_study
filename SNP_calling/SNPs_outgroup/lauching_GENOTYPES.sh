#!/bin/bash
#SBATCH -J conc
#SBATCH -o conc.%j.out
#SBATCH -e conc.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH -c 1
#SBATCH --mem 200G

module load python/3.9

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr1.txt -w GENOTYPES_SRR7440095_chr1.txt -o GENOTYPES_SRR7440094_95_CONCATENATED_chr1.txt










