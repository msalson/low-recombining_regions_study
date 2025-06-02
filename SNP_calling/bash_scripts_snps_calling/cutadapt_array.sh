#!/bin/bash
#SBATCH -J cutadapt
#SBATCH -e array_%A-%a.err
#SBATCH -A snp_calling_gatk_mil      
#SBATCH --mail-user marine.salson@ird.fr
#SBATCH --mail-type=FAIl,END
#SBATCH -p fast
#SBATCH --array=1-15%2
#SBATCH --cpus-per-task=10          # Allocation de 1 CPU par Task
#SBATCH --output=array_%A-%a.out

module load cutadapt/3.1

#to get the sample name from text file:
indiv_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ID_simple_index.txt) 

#{indiv_id}=IP8767
#fastq_data/${indiv_id}*1.fastq.gz
#fastq_data/${indiv_id}*2.fastq.gz

cutadapt -j 10 -m 35 -q 20 -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -b agatcggaagagcggttcagcaggaatgccgagaccgatctcgtatgccgtcttctgcttg -b caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct -b tctagccttctcgccaagtcgtccttacggctctggctagagcatacggcagaagacgaac -b gttcgtcttctgccgtatgctctagccagagccgtaaggacgacttggcgagaaggctaga -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -b CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -b TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTGNNNNNNTAGAGCATACGGCAGAAGACGAAC -b GTTCGTCTTCTGCCGTATGCTCTANNNNNNCACTGACCTCAAGTCTGCACACGAGAAGGCTAGA -b AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTC -b GAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT -b CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTgctcttccgatct -b agatcggaagagcACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG -b ttactatgccgctggtggctctagatgtgagaaagggatgtgctgcgagaaggctaga -b tctagccttctcgcagcacatccctttctcacatctagagccaccagcggcatagtaa -B AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -B agatcggaagagcggttcagcaggaatgccgagaccgatctcgtatgccgtcttctgcttg -B caagcagaagacggcatacgagatcggtctcggcattcctgctgaaccgctcttccgatct -B tctagccttctcgccaagtcgtccttacggctctggctagagcatacggcagaagacgaac -B gttcgtcttctgccgtatgctctagccagagccgtaaggacgacttggcgagaaggctaga -B AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -B CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -B TCTAGCCTTCTCGTGTGCAGACTTGAGGTCAGTGNNNNNNTAGAGCATACGGCAGAAGACGAAC -B GTTCGTCTTCTGCCGTATGCTCTANNNNNNCACTGACCTCAAGTCTGCACACGAGAAGGCTAGA -B AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTC -B GAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT -B CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTgctcttccgatct -B agatcggaagagcACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG -B ttactatgccgctggtggctctagatgtgagaaagggatgtgctgcgagaaggctaga -B tctagccttctcgcagcacatccctttctcacatctagagccaccagcggcatagtaa -o ${indiv_id}_CU_1.fastq.gz -p ${indiv_id}_CU_2.fastq.gz ${indiv_id}*1.fastq.gz ${indiv_id}*2.fastq.gz 2> cutadapt_log/${indiv_id}.err > cutadapt_log/${indiv_id}.log
