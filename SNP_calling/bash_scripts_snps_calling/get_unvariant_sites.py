#!/usr/bin/env python3

import os,re
import argparse

parser = argparse.ArgumentParser(description='Script to get invariant sites')
parser.add_argument('-v', '--vcf_file', help='vcf file', required=True)

args = parser.parse_args()

vcf=args.vcf_file

with open("chr1_non_variant_sites.vcf", "w") as fR:
	with open (vcf, "r", encoding="utf8", errors='ignore') as file:
		for line in file:
			e=line.split("\t")
			if str(e[0][0]) != "#":
				alt=e[4]
				if str(alt) == ".": 
					c=0 
					for i in range(9,135):
						search=e[i].split(":") #GT:DP:RGQ
						if len(search)>1:
							GT=search[0]
							DP=search[1]
							if int(DP)>=5 and int(DP)<100 and str(GT) == "0/0" :
								c+=1
					if int(c) >= 63:
						fR.write(line)


