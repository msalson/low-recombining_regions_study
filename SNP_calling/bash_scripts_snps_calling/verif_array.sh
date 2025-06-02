#!/bin/bash
#SBATCH -J verif_Super-Scaffold_XX
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p long
#SBATCH --array=1-178%2
#SBATCH --cpus-per-task=1          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out

# To check if the gvcf has been written until the end

indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../ID_all.txt) 

zcat ${indiv_id}_Super-Scaffold_XX.gvcf.gz | tail -n 1 > verif/${indiv_id}.txt

