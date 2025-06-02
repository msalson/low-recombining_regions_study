#!/bin/bash
#SBATCH -J bwa_array
#SBATCH -o bwa."%j"-"%a".out
#SBATCH -e bwa."%j"-"%a".err
#SBATCH --mail-user marine.salson@etu.umontpellier.fr
#SBATCH --mail-type=ALL
#SBATCH -p highmem
#SBATCH -w node7
#SBATCH -c 40
#SBATCH --array=1-130%10
#SBATCH --cpus-per-task=10          

module load bwa-mem2/2.2.1
module load picard/2.5.0 
module load samtools/1.9

indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.txt) 

#${indiv_id}=trimfiltcutI1I8-T7
#trimfiltcutI10I11-T37_R1_paired.fastq => ${indiv_id}_R1_paired.fastq
#trimfiltcutI10I11-T37_R2_paired.fastq => ${indiv_id}_R2_paired.fastq

#bwa mem :

#bwa-mem2 index reference.fasta

bwa-mem2 mem -t 10 reference.fasta ${indiv_id}_R1_paired.fastq ${indiv_id}_R2_paired.fastq -R '@RG\tID:'+${indiv_id}+'\tSM:'+${indiv_id}+'\tPL:Illumina' -o results_mapping/${indiv_id}.sam

#picard tools:
java -jar /usr/local/bioinfo/picard-tools-2.5.0/picard.jar SortSam SORT_ORDER=coordinate TMP_DIR=tmp VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE  INPUT=results_mapping/${indiv_id}.sam OUTPUT=results_mapping/${indiv_id}.PICARDTOOLSSORT.bam

java -jar /usr/local/bioinfo/picard-tools-2.5.0/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=results_mapping/${indiv_id}.PICARDTOOLSSORT.bam O=results_mapping/${indiv_id}.PICARDTOOLSSORT.dedup.bam M=results_mapping/${indiv_id}_metrices_MK.txt

#samtools

samtools flagstat results_mapping/${indiv_id}.PICARDTOOLSSORT.dedup.bam > samtools/${indiv_id}_samtoolsflagstat.log 2> samtools/${indiv_id}_samtoolsflagstat.err

samtools view -b -h -@ 40 -f 0x02 -o results_mapping/${indiv_id}.SAMTOOLSVIEW.bam results_mapping/${indiv_id}.PICARDTOOLSSORT.dedup.bam

samtools index results_mapping/${indiv_id}.SAMTOOLSVIEW.bam

