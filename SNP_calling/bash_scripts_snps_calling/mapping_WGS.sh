#!/bin/bash
#SBATCH -J mapping
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH --array=1-173%10
#SBATCH --cpus-per-task=50          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out

module load bwa-mem2/2.2.1
module load picard/2.18.9 
module load samtools/1.9

#obtention of sample name in text file:
indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ID_all.txt) 

#${indiv_id}=IP8767
#${indiv_id}*1.fastq.gz
#${indiv_id}*2.fastq.gz

#bwa mem :
bwa-mem2 mem -t 50 pearl_millet_ONT.fasta ./fastq/${indiv_id}*_1.fastq.gz ./${indiv_id}*_2.fastq.gz -R '@RG\tID:${indiv_id}\tSM:${indiv_id}\tPL:Illumina' -o ./results_bwa/${indiv_id}.sam

#picard tools sortsam :
picard SortSam SORT_ORDER=coordinate TMP_DIR=tmp VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE INPUT=./results_bwa/${indiv_id}.sam OUTPUT=./results_bwa/${indiv_id}.PICARDTOOLSSORT.bam

samtools flagstat ./results_bwa/${indiv_id}.PICARDTOOLSSORT.bam > ./resultats_flagstat/${indiv_id}_samtoolsflagstat.log 2> ./resultats_flagstat/${indiv_id}_samtoolsflagstat.err

rm ./results_bwa/${indiv_id}.sam

#samtools

samtools view -b -h -@ 50 -f 0x02 -o ./results_bwa/${indiv_id}.SAMTOOLSVIEW.bam ./results_bwa/${indiv_id}.PICARDTOOLSSORT.bam

rm ./results_bwa/${indiv_id}.PICARDTOOLSSORT.bam

samtools flagstat ./results_bwa/${indiv_id}.SAMTOOLSVIEW.bam > ./resultats_flagstat/${indiv_id}_VIEW_samtoolsflagstat.log 2> ./resultats_flagstat/${indiv_id}_VIEW_samtoolsflagstat.err

