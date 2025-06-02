#!/bin/bash
#SBATCH -J outgroup
#SBATCH -o outgroup.%j.out
#SBATCH -e outgroup.%j.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH -c 1
#SBATCH --mem 80G

module load python/3.9

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr1_full_all_reads.txt -w GENOTYPES_SRR7440095_chr1_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr1_full.txt

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr2_full_all_reads.txt -w GENOTYPES_SRR7440095_chr2_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr2_full.txt

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr4_full_all_reads.txt -w GENOTYPES_SRR7440095_chr4_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr4_full.txt

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr5_full_all_reads.txt -w GENOTYPES_SRR7440095_chr5_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr5_full.txt

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr6_full_all_reads.txt -w GENOTYPES_SRR7440095_chr6_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr6_full.txt

./GENOTYPES_outgroup_concatenate.py -v GENOTYPES_SRR7440094_chr7_full_all_reads.txt -w GENOTYPES_SRR7440095_chr7_full_all_reads.txt -o GENOTYPES_SRR7440094-5_CONCATENATED_chr7_full.txt

