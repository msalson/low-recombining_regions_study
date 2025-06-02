#!/usr/bin/env python3

import os,re
import argparse

parser = argparse.ArgumentParser(description='Script to get genomic windows positions.')
parser.add_argument('-v', '--vcf_file', help='vcf file', required=True)
parser.add_argument('-w', '--windows_size', help='number of snps per windows', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)

args = parser.parse_args()

VCF = args.vcf_file
windowsSize = int(args.windows_size)
outputFile = args.output

c=0
formerChromosome=""
list_chr=[]
dico_chr={}

with open(VCF, "r", encoding="utf8", errors='ignore') as vcfFile :
    for line in vcfFile :
        if line[0]!="#":
            e=line.split("\t")
            chromosome=e[0]
            if str(chromosome) not in list_chr:
                list_chr.append(chromosome)
            if str(chromosome) not in dico_chr:   
                dico_chr[chromosome]={}
                dico_chr[chromosome]["posB"]=[]
                dico_chr[chromosome]["posE"]=[]
            if str(chromosome)!=str(formerChromosome) :
                c=0
            c+=1
            position=e[1]
            if int(c%windowsSize) == int(1) :
                posStart=position
            if int((c+1)%windowsSize) == int(1) :
                posEnd=position
                dico_chr[chromosome]["posB"].append(posStart)
                dico_chr[chromosome]["posE"].append(posEnd)
            formerChromosome=e[0]
                    
sorted_list_chr = sorted(list_chr) # alphabetic order as lostruct handles chromosomes in vcf file 

id=0
with open(outputFile, "w") as Rfile :
    Rfile.write("ID\tChr\tposB\tposE\n")
    for chromosome in sorted_list_chr: #['chr1', 'chr10', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']
        nb_pos=len(dico_chr[chromosome]["posB"])
        for i in range(0,int(nb_pos)):
            id+=1
            posStart=dico_chr[chromosome]["posB"][i]
            posEnd=dico_chr[chromosome]["posE"][i]
            Rfile.write(str(id)+"\t"+str(chromosome)+"\t"+str(posStart)+"\t"+str(posEnd)+"\n")

