#!/usr/bin/env python3

import os,re
import argparse
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser(description='Script to generate a table of genotypes.')
parser.add_argument('-v', '--vcf_file', help='VCF file with the samples ordered', required=True)
parser.add_argument('-s', '--outlier_snps', help='output of the PCAdapt analysis, with the outlier SNPs', required=True)
parser.add_argument('-m', '--output_matrix', help='matrix of genotypes', required=True)
parser.add_argument('-p', '--output_pos', help='file with positions of the outlier SNPs for x axis of the table of genotype', required=True)
parser.add_argument('-f', '--output_png_file', help='final png file with the table of genotypes', required=True)

args = parser.parse_args()

VCF = args.vcf_file

vcf_ordered=args.vcf_file
file_with_outlier_SNPs=args.outlier_snps 
output_matrix_geno=args.output_matrix
output_index_pos_SNPs_outliers=args.output_pos
output_png_file=args.output_png_file

# To get the positions of outlier SNPs identified with PCAdapt:
list_outlier_SNPs=[]
with open(file_with_outlier_SNPs, "r") as fileSNPs :
    for line in fileSNPs:
        e=line.split(" ")
        outlier=e[1]
        list_outlier_SNPs.append(outlier)

# To get the matrix of genotypes:
count_snp=0
c=0
with open(output_matrix_geno, "w") as fmatrix:
    with open(vcf_ordered, "r", encoding="utf8", errors='ignore') as vcfFile :
        for line in vcfFile :
            e=line.split("\t")
            if str(e[0]) == "#CHROM": 
                nb_columns=len(e)
                nb_samples=int(nb_columns)-9
            if line[0]!="#":
                count_snp+=1
                if str(count_snp) in list_outlier_SNPs:
                    c+=1
                    for i in range(9,int(nb_columns)):
                        infos=e[i]
                        genotypeSearch=infos.split(":")
                        geno=genotypeSearch[0]
                        if str(geno)=="./.":
                            geno=1.0
                        if str(geno)=="0/1" or str(geno)=="0|1":
                            geno=3.0
                        if str(geno)=="1/1" or str(geno)=="1|1":
                            geno=2.0
                        if str(geno)=="0/0" or str(geno)=="0|0":
                            geno=4.0
                        if int(i) < (int(nb_columns)-1):
                            fmatrix.write(str(geno)+"\t")
                        if int(i) == (int(nb_columns)-1):
                            fmatrix.write(str(geno))
                fmatrix.write("\n")

data = pd.read_csv(output_matrix_geno, sep="\t", header=None)
data = data.transpose()


# To get the positions of SNPs for the x axis of the matrix:
i=0
with open(output_index_pos_SNPs_outliers, "w") as fR:
    fR.write("chr\tnb\tpos\n")
    with open(vcf_ordered, "r", encoding="utf8", errors='ignore') as vcfFile :
        for line in vcfFile :
            e=line.split("\t")
            if str(e[0][0]) != "#": 
                chromosome=e[1]
                position=e[2]
                i+=1
                if str(i) in list_outlier_SNPs:
                    fR.write(str(chromosome)+"\t"+str(c)+"\t"+str(position)+"\n")


list_pos_SNPs_matrix_sep=[]
c=0
nb_snps=1000
with open(output_index_pos_SNPs_outliers, "r") as f:
    for line in f:
        e=line.split("\t")
        if str(e[0])!="chr":
            c+=1
            position_SNP=e[0]
            new_pos=round((int(position_SNP)/1000000),1)
            if c==1:
                list_pos_SNPs_matrix_sep.append(new_pos)
            if c==int(nb_snps): 
                list_pos_SNPs_matrix_sep.append(new_pos) 
            if c%15 == 0:
                list_pos_SNPs_matrix_sep.append(new_pos)
            else:
                if c!=1 and c!= int(nb_snps):
                    list_pos_SNPs_matrix_sep.append("")

x_axis_labels = list_pos_SNPs_matrix_sep

# Obtention of the heatmap:
myColors = ((0.0, 0.0, 0.0, 0.1), (0.8, 0.0, 0.0, 1.0), (0.3, 0.0, 0.8, 0.30), (0.0, 0.2, 0.8, 0.8))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
heat_map = sb.heatmap(data,cmap=cmap,yticklabels=False,xticklabels=x_axis_labels) #,cmap="Blues"
heat_map.tick_params(left=False, bottom=False)
colorbar = heat_map.collections[0].colorbar
colorbar.set_ticks([1.4, 2.1, 2.9, 3.6])
colorbar.set_ticklabels(['missing data', ' 1 / 1', ' 0 / 1',' 0 / 0'])
heat_map.set_ylabel(str(nb_samples)+' samples')
heat_map.set_xlabel('SNPs')
heat_map.set_xticklabels(x_axis_labels, rotation = 35, fontsize = 7)

heat_map.figure.savefig(output_png_file,dpi=600)
