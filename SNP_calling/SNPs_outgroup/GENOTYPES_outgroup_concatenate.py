#!/usr/bin/env python3

import os,re
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Script to get genotypes of the outgroups of pearl millet from gvcf file (generated with BP_RESOLUTION parameter).')
parser.add_argument('-v', '--file1', help='Gvcf file', required=True)
parser.add_argument('-w', '--file2', help='Gvcf file', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)

args = parser.parse_args()

file94 = args.file1
file95 = args.file2
outputFile = args.output

dico_94={}
dico_95={}

with open(file94, "r") as f94:
    for line in f94:
        line_split=line.split("\t")
        position=line_split[1]
        if str(position)!=str("positions"):
            ref=line_split[2]
            sampleGT=line_split[3]
            dp=line_split[4]
            gq=line_split[5][:-1]
            dico_94[position]={}
            dico_94[position]["sampleGT94"]=sampleGT
            dico_94[position]["DP94"]=dp
            dico_94[position]["GQ94"]=gq
            dico_94[position]["REF"]=ref

with open(file95, "r") as f95:
    for line in f95:
        line_split=line.split("\t")
        position=line_split[1]
        if str(position)!=str("positions"):
            ref=line_split[2]
            sampleGT=line_split[3]
            dp=line_split[4]
            gq=line_split[5][:-1]
            dico_95[position]={}
            dico_95[position]["sampleGT95"]=sampleGT
            dico_95[position]["DP95"]=dp
            dico_95[position]["GQ95"]=gq
            dico_95[position]["REF"]=ref

with open(outputFile, "w") as file_R:
    file_R.write("position\tREF\tGT1\t\tDP1\tGQ1\tGT2\tDP2\tGQ2\n")
    for i in range(1,315000000):
        if str(i) in dico_94:
            position=i
            REF=dico_94[str(i)]["REF"]
            gt94=dico_94[str(i)]["sampleGT94"]
            dp94=dico_94[str(i)]["DP94"]
            gq94=dico_94[str(i)]["GQ94"]
        else :
            gt94="."
            dp94="."
            gq94="."     
        if str(i) in dico_95:
            position=i
            REF=dico_95[str(i)]["REF"]
            gt95=dico_95[str(i)]["sampleGT95"]
            dp95=dico_95[str(i)]["DP95"]
            gq95=dico_95[str(i)]["GQ95"]
        else:
            gt95="."
            dp95="."
            gq95="."    
        if (gt95!=".") or (gt94!="."):
            file_R.write(str(position)+"\t"+str(REF)+"\t"+str(gt94)+"\t"+str(dp94)+"\t"+str(gq94)+"\t"+str(gt95)+"\t"+str(dp95)+"\t"+str(gq95)+"\n")
