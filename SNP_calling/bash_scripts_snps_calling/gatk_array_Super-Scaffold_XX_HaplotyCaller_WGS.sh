#!/bin/bash
#SBATCH -J gatk_Super-Scaffold_XX
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH --array=1-178%45
#SBATCH --cpus-per-task=1          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out
#SBATCH --mem-per-cpu 15G

module load gatk4/4.2.3.0

indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./ID_all.txt)

gatk HaplotypeCaller --java-options -Xmx14G -R pearl_millet_ONT.fasta -I ./results_bwa/${indiv_id}.bam -O ./Super-Scaffold_XX/${indiv_id}_Super-Scaffold_XX.gvcf.gz -L Super-Scaffold_XX_polished --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF

