#!/usr/bin/env python3

# This script was made to transfer the positions of the SNPs from the scaffolds to the chromosomes of the new pearl millet genome (Salson et al. 2023).
# A file with the orientation and the start and end positions of the scaffolds within each chromosome was used.

# We checked if the positions of the SNPs on the chromosomes were correct:
# We checked if the REF alleles of the SNPs on the new VCF file were the same as the nucleotides on the sequences of the chromosome for the same positions. 
# We also used the 76,018 poisitions on the vcf file of the spatial dataset for which the mapping and SNPs calling were made directly to the chromosomes of the new pearl millet genome;
# for each of the common snps, we checked if the ALT alleles were the same: it was 100% correct. 


import os,re
import argparse

parser = argparse.ArgumentParser(description='Script to transfer the SNPs of a VCF file from the scaffolds to chromosomes.')
parser.add_argument('-v', '--vcf_file', help='vcf file with snps on scaffolds', required=True)
parser.add_argument('-w', '--file_infos', help='file with start and end positions of the scaffolds within each chromosome', required=True)
parser.add_argument('-o', '--output', help='output vcf file with snps on chromosomes', required=True)

args = parser.parse_args()

VCF = args.vcf_file
infos_file = args.file_infos
outputFile = args.output

dico={}
with open(infos_file, "r") as csv:
    for line in csv:
        i=line.split("\t")
        scaffold=i[0]
        chromosome=i[1]
        orientation=i[4]
        length=i[5][:-1]
        start=i[2]
        end=i[3]
        if str(scaffold)!="scaffolds":
            dico[scaffold]={}
            dico[scaffold]["chr"]=chromosome
            dico[scaffold]["orientation"]=orientation
            dico[scaffold]["start"]=start
            dico[scaffold]["end"]=end
            dico[scaffold]["length"]=length

scaffolds_on_chr=[]
for scaffold in dico:
    scaffolds_on_chr.append(scaffold)

with open(outputFile, "w") as fR:
    with open(VCF, "r", encoding="utf8", errors='ignore') as vcfFile :
        for line in vcfFile :
            e=line.split("\t")
            if str(e[0][0]) != "#": 
                scaffold=e[0]
                if str(scaffold) in scaffolds_on_chr:
                    position_SNP_on_the_scaffold=e[1]
                    o=dico[scaffold]["orientation"]
                    allele_ref=e[3]
                    allele_ref_new=e[3]
                    new_chr="chr"+str(dico[scaffold]["chr"])
                    position_start_scaffold=dico[scaffold]["start"]
                    new_pos=int(position_start_scaffold)+int(position_SNP_on_the_scaffold)-1
                    if str(o) == "-": 
                        if str(allele_ref)=="A":
                            allele_ref_new="T"
                        if str(allele_ref)=="T":
                            allele_ref_new="A"
                        if str(allele_ref)=="G":
                            allele_ref_new="C"
                        if str(allele_ref)=="C":
                            allele_ref_new="G"
                        length_saffold=dico[scaffold]["length"]
                        new_pos_SNP=int(length_saffold)-int(position_SNP_on_the_scaffold)
                        new_pos=int(position_start_scaffold)+(new_pos_SNP)
                    allele_alt=e[4]
                    allele_alt_final=e[4]
                    if str(o) == "-": 
                        if str(allele_alt)=="A":
                            allele_alt_final="T"
                        if str(allele_alt)=="T":
                            allele_alt_final="A"
                        if str(allele_alt)=="G":
                            allele_alt_final="C"
                        if str(allele_alt)=="C":
                            allele_alt_final="G"
                    for i in range(0,182):
                        if str(o) == "-": 
                            if i<181 and (i!=0) and (i!=1) and (i!=2) and (i!=3) and (i!=4) :
                                fR.write(str(e[i])+"\t")
                            if i==181 :
                                fR.write(str(e[i]))
                            if i == 0:
                                fR.write(str(new_chr)+"\t")
                            if i == 1:
                                fR.write(str(new_pos)+"\t")
                            if i == 2:
                                fR.write(".\t"+str(allele_ref_new)+"\t")
                            if i == 3:
                                fR.write(str(allele_alt_final)+"\t")
                        if str(o) == "+": 
                            if i<181 and (i!=0) and (i!=1) :
                                fR.write(str(e[i])+"\t")
                            if i==181 and (i!=0) and (i!=1) :
                                fR.write(str(e[i]))
                            if i == 0:
                                fR.write(str(new_chr)+"\t")
                            if i == 1:
                                fR.write(str(new_pos)+"\t")
                    if str(scaffold) not in scaffolds_on_chr:
                        for i in range(0,182):
                            if i<181 :
                                fR.write(str(e[i])+"\t")
                            if i==181 :
                                fR.write(str(e[i]))

