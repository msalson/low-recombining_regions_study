#!/bin/bash
#SBATCH -J gatk_SRR7440095
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH --array=1-1%1
#SBATCH --cpus-per-task=10          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out

module load gatk4/4.2.3.0
module load samtools/1.9

gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr1_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr1 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION
#gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr2_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr2 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION
#gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr4_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr4 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION
#gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr5_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr5 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION
#gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr6_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr6 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION
#gatk HaplotypeCaller --java-options -Xmx30G -R pearl_millet_23DB_ONT_assembly_all_chr_and_unplaced_sequences.fasta -I SRR7440095_PICARDTOOLSSORT_f.bam -O SRR7440095_chr7_BP_RESOLUTION_f_all_reads.gvcf.gz -L chr7 --output-mode EMIT_ALL_ACTIVE_SITES -ERC BP_RESOLUTION

module load python/3.9

gunzip SRR7440095_chr1_BP_RESOLUTION_f_all_reads.gvcf.gz
#gunzip SRR7440095_chr2_BP_RESOLUTION_f_all_reads.gvcf.gz
#gunzip SRR7440095_chr4_BP_RESOLUTION_f_all_reads.gvcf.gz
#gunzip SRR7440095_chr5_BP_RESOLUTION_f_all_reads.gvcf.gz
#gunzip SRR7440095_chr6_BP_RESOLUTION_f_all_reads.gvcf.gz
#gunzip SRR7440095_chr7_BP_RESOLUTION_f_all_reads.gvcf.gz

./GENOTYPES_outgroup_chr.py -v SRR7440095_chr1_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr1_full_all_reads.txt
#./GENOTYPES_outgroup_chr.py -v SRR7440095_chr2_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr2_full_all_reads.txt
#./GENOTYPES_outgroup_chr.py -v SRR7440095_chr4_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr4_full_all_reads.txt
#./GENOTYPES_outgroup_chr.py -v SRR7440095_chr5_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr5_full_all_reads.txt
#./GENOTYPES_outgroup_chr.py -v SRR7440095_chr6_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr6_full_all_reads.txt
#./GENOTYPES_outgroup_chr.py -v SRR7440095_chr7_BP_RESOLUTION_f_all_reads.gvcf -o results_Maud/GENOTYPES_SRR7440095_chr7_full_all_reads.txt

gzip SRR7440095_chr1_BP_RESOLUTION_f_all_reads.gvcf.gz
#gzip SRR7440095_chr2_BP_RESOLUTION_f_all_reads.gvcf.gz
#gzip SRR7440095_chr4_BP_RESOLUTION_f_all_reads.gvcf.gz
#gzip SRR7440095_chr5_BP_RESOLUTION_f_all_reads.gvcf.gz
#gzip SRR7440095_chr6_BP_RESOLUTION_f_all_reads.gvcf.gz
#gzip SRR7440095_chr7_BP_RESOLUTION_f_all_reads.gvcf.gz



