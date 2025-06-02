#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
import os
import argparse


parser = argparse.ArgumentParser(description='Script to get diploid fasta sequences from GFF and ALT sequences')

#previous step:
#module load bcftools/1.9
#module load gatk4/4.2.3.0
#bcftools query -l vcf_126_mils_76018_filtered_SNPs_annotated.vcf  | while read S
#do
#gatk FastaAlternateReferenceMaker -R pearl_millet_23DB_ONT_assembly.fasta -O "${S}_ALL_CHR_FULL_IUPAC.fasta" --use-iupac-sample "${S}" -V vcf_126_mils_76018_filtered_SNPs_annotated.vcf
#done

parser.add_argument('-g', '--gff_cds', help='gff with CDS', required=True)
parser.add_argument('-i', '--indiv', help='individual identifier', required=True)
args = parser.parse_args()

GFF= args.gff_cds
individual=args.indiv

# So-21-30218-02_ALL_CHR_FULL_IUPAC.fasta
# 126 files as:
# *_ALL_CHR_FULL_IUPAC.fasta

####################################
## Get dictionary with ALT sequences
####################################

file_ALT = str(individual)+"_ALL_CHR_FULL_IUPAC.fasta"

new_file_ALT = str(individual)+"_ALL_CHR_FULL_IUPAC_HEADER.fasta"

with open(new_file_ALT, "w") as new_fasta:
    for record in SeqIO.parse(file_ALT, "fasta"):
        header=str(">chr"+str(record.id)+"\n")
        new_fasta.write(header)
        seq=str(record.seq)
        new_fasta.write(seq+"\n")

dico_chr_ALT={}

for record in SeqIO.parse(new_file_ALT, "fasta"):
    if str(record.id) not in dico_chr_ALT:
        dico_chr_ALT[str(record.id)]=record.seq

##########################################################
## Extract CDS from GFF file considering phase information
##########################################################

fasta_file_REVERSE=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_reverse.fasta"
fasta_file_FORWARD=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_forward.fasta"

count=1
count_reverse=0
list_CDS=[]
with open(fasta_file_REVERSE, "w") as fasta_to_generate_rev:
    with open(fasta_file_FORWARD, "w") as fasta_to_generate_for:
        with open(GFF, 'r') as gff:
            for line in gff:
                e=line.split("\t")
                chromosome=e[0]
                forward_reverse=e[6]
                if str(chromosome) in dico_chr_ALT: 
                    start=e[3] 
                    end=e[4] 
                    phase=e[7]
                    sequence=dico_chr_ALT[str(chromosome)][int(start)-1:int(end)]  
                    check_doublon=str(sequence)+str(start)+str(end)
                    if str(check_doublon) not in list_CDS:
                        list_CDS.append(check_doublon)
                        
                        if str(forward_reverse) == "-":    
                            if str(phase) == "0" : 
                                sequence=sequence.reverse_complement()
                                fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                            elif str(phase) == "1" :
                                sequence=sequence.reverse_complement()
                                sequence=sequence[1:]
                                fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                            elif str(phase) == "2" : 
                                sequence=sequence.reverse_complement()
                                sequence=sequence[2:]
                                fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                                                          
                        if str(forward_reverse) == "+":   
                            if str(phase) == "0" : 
                                fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                            elif str(phase) == "1" : 
                                sequence=sequence[1:]
                                fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                            elif str(phase) == "2" : 
                                sequence=sequence[2:]
                                fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                                                                         
                    count+=1

#############################################################
## Obtain diploid sequences for each individual and each CDS
#############################################################

# Forward CDS:

file_to_duplicate_seq_for=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_forward.fasta"
file_duplicated_forward=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_forward_DIPLOID.fasta"

no_change=['A','T','G','C']

count=1
with open(file_duplicated_forward, "w") as fR:
    for record in SeqIO.parse(str(file_to_duplicate_seq_for), "fasta"):
        sequence_1=""
        sequence_2=""
        id_seq=str(record.id) #>CDS85_chr1_3621451_3621710
        id_seq_e=id_seq.split("_")
        header=str(id_seq_e[1])+"_"+str(id_seq_e[2])+"_"+str(id_seq_e[3])
        SEQ=record.seq
        for car in SEQ:
            if str(car) == ".":
                sequence_1=str(sequence_1)+"N"
                sequence_2=str(sequence_2)+"N"
            if str(car) in no_change:
                sequence_1=str(sequence_1)+str(car)
                sequence_2=str(sequence_2)+str(car)
            if str(car) == "W":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "S":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "M":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"C"
            if str(car) == "K":
                sequence_1=str(sequence_1)+"G"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "R":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "Y":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"T"
        if  len(sequence_1) == len(sequence_2) and len(sequence_1) == len(SEQ): #check that sequences have the same length, as expected
            #>CDS_1|Pearl_millet|So-21-28371-02|Allele_1
            header_1=">CDS_"+str(count)+"_"+str(header)+"|Pearl_millet|"+str(individual)+"|Allele_1"
            fR.write(str(header_1)+"\n")
            fR.write(str(sequence_1)+"\n")
            header_2=">CDS_"+str(count)+"_"+str(header)+"|Pearl_millet|"+str(individual)+"|Allele_2"
            fR.write(str(header_2)+"\n")
            fR.write((sequence_2)+"\n")
        else:
            print(header)
        count+=1

# Reverse CDS:

file_to_duplicate_seq_rev=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_MANUAL_reverse.fasta"
file_duplicated_reverse=str(individual)+"_ALL_CHR_FULL_IUPAC_CDS_MANUAL_reverse_DIPLOID.fasta"

count=1
with open(file_duplicated_reverse, "w") as fR:
    for record in SeqIO.parse(str(file_to_duplicate_seq_rev), "fasta"):
        sequence_1=""
        sequence_2=""
        id_seq=str(record.id) #>CDS85_chr1_3621451_3621710
        id_seq_e=id_seq.split("_")
        header=str(id_seq_e[1])+"_"+str(id_seq_e[2])+"_"+str(id_seq_e[3])
        SEQ=record.seq
        for car in SEQ:
            if str(car) == ".":
                sequence_1=str(sequence_1)+"N"
                sequence_2=str(sequence_2)+"N"
            if str(car) in no_change:
                sequence_1=str(sequence_1)+str(car)
                sequence_2=str(sequence_2)+str(car)
            if str(car) == "W":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "S":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "M":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"C"
            if str(car) == "K":
                sequence_1=str(sequence_1)+"G"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "R":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "Y":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"T"
        if  len(sequence_1) == len(sequence_2) and len(sequence_1) == len(SEQ): #check that sequences have the same length, as expected
            #>CDS_1|Pearl_millet|So-21-28371-02|Allele_1
            header_1=">CDS_"+str(count)+"_"+str(header)+"|Pearl_millet|"+str(individual)+"|Allele_1"
            fR.write(str(header_1)+"\n")
            fR.write(str(sequence_1)+"\n")
            header_2=">CDS_"+str(count)+"_"+str(header)+"|Pearl_millet|"+str(individual)+"|Allele_2"
            fR.write(str(header_2)+"\n")
            fR.write((sequence_2)+"\n")
        else:
            print(header)
        count+=1

os.system('mv '+str(fasta_file_REVERSE)+' results/')
os.system('mv '+str(fasta_file_FORWARD)+' results/')
os.system('mv '+str(file_duplicated_forward)+' results/')
os.system('mv '+str(file_duplicated_reverse)+' results/')


