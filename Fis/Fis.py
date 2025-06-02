#!/usr/bin/env python3

import os,re
import argparse

# The script required access to VCFtools. The version 0.1.17 was used in the study. 

parser = argparse.ArgumentParser(description='Script to get Fis for each locus from a VCF file.')
parser.add_argument('-v', '--vcf_file', help='vcf file', required=True)
parser.add_argument('-o', '--output_file', help='output file name', required=True)

args = parser.parse_args()

VCF = args.vcf_file
output = args.output_file

suffix=str(VCF).split(".")
suffix=suffix[len(suffix)-1]

if str(suffix) == "vcf":
	os.system('vcftools --vcf '+str(VCF)+' --hardy')
if str(suffix) == "gz":
	os.system('vcftools --gzvcf '+str(VCF)+' --hardy')

os.system('mv out.hwe '+str(VCF)+'.hwe')

file_hwe=str(VCF)+".hwe" #ex in the study: spatial_dataset_126_samples.vcf.hwe, accessible in the repository. 

with open(output, "w") as Fr:
	Fr.write("chromosome\tlocus\tFis\theteroObs\theteroExp\n")
	with open(file_hwe, "r") as f:
		for line in f:
			e=line.split("\t")
			if str(e[0])!="CHR":
				chromosome=e[0] #chr1
				locus=e[1] #36707
				infos1=e[2] #24/24/21 => OBS(HOM1/HET/HOM2)
				infos2=e[3] #18.78/34.43/15.78 => E(HOM1/HET/HOM2)
				i1=infos1.split("/")            
				i2=infos2.split("/")
				heteroO=float(i1[1]) #24
				heteroE=float(i2[1]) #34.43
				if int(heteroE)>0:
					Fis=(heteroE-heteroO)/heteroE
					Fr.write(str(chromosome)+"\t"+str(locus)+"\t"+str(Fis)+"\t"+str(heteroO)+"\t"+str(heteroE)+"\n")
