#!/bin/bash
#SBATCH -J gatk
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH --array=1-126%20
#SBATCH --cpus-per-task=5          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out
#SBATCH --mem-per-cpu 10G

module load gatk4/4.2.3.0

indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./samples.txt)

gatk HaplotypeCaller --java-options -Xmx9G -R reference.fasta -I ${indiv_id}.MK.SAMTOOLSVIEW.bam -O chr1_output/${indiv_id}_chr1.gvcf.gz -L chr1 --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF

