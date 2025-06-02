#!/usr/bin/env python3

import os,re
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Script to get genotypes of the outgroups of pearl millet from gvcf file (generated with BP_RESOLUTION parameter).')
parser.add_argument('-v', '--gvcf_file', help='Gvcf file', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)

args = parser.parse_args()

file_vcf = args.gvcf_file
outputFile = args.output


sample_GT=""
with open(outputFile, "w") as fR:
    fR.write("chromosome\tpositions\tREF\tsample_GT\tDP\tGQ\n")
    with open(file_vcf, "r", encoding="utf8", errors='ignore') as vcfFile :
        for line in vcfFile :
            e=line.split("\t")
            if str(e[0][0]) != "#": 
                new_chr=e[0]
                new_pos=e[1]
                allele_ref_new=e[3]
                if (str(allele_ref_new)=="A" or str(allele_ref_new)=="T" or str(allele_ref_new)=="G" or str(allele_ref_new)=="C"):
                    allele_alt=e[4] #G,<NON_REF>
                    if str(allele_alt) == "<NON_REF>":
                        info_field=e[8] #100% GT:AD:DP:GQ:PL
                        info_field_sample=e[9]#0/0:1:3:1:0,3,35
                        info_field_sample_split=info_field_sample.split(":")
                        GT=info_field_sample_split[0]
                        GQ=info_field_sample_split[3]
                        DP=info_field_sample_split[2]
                        if int(GQ)>=1 and int(DP)>=1:
                            sample_GT=str(allele_ref_new)+str(allele_ref_new)
                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                    if str(allele_alt) != "<NON_REF>" :
                        allele_alt_split=allele_alt.split(",") #G,<NON_REF>
                    ####################                   
                    #####G,<NON_REF>
                    ####################
                        if len(allele_alt_split) == 2:#G,<NON_REF>
                            ALT1=allele_alt_split[0]
                            info_field=e[8] #4 possibilite GT:AD:DP:GQ:PL, etc
                            info_field_split=info_field.split(":") #str(info_field_split[1])=="AD"
                            if str(info_field_split[0])=="GT" and str(info_field_split[1])=="AD" and str(info_field_split[2])=="DP" and str(info_field_split[3])=="GQ":
                                info_field_sample=e[9]#0/0:1:3:1:0,3,35
                                info_field_sample_split=info_field_sample.split(":")
                                GT=info_field_sample_split[0]
                                GQ=info_field_sample_split[3]
                                DP=info_field_sample_split[2]
                                if (str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C") :
                                    if (str(GT)=="0/0" or str(GT)=="0|0") and int(GQ)>=1 and int(DP)>=1:
                                        sample_GT=str(allele_ref_new)+str(allele_ref_new)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                    if (str(GT)=="1/1" or str(GT)=="1|1") and int(GQ)>=1 and int(DP)>=1:
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(new_ALT1)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                    if (str(GT)=="0/1" or str(GT)=="0|1" or str(GT)=="1|0") and int(GQ)>=1 and int(DP)>=1:
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(allele_ref_new)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                    ####################
                    #####A, G,<NON_REF>                                    
                    ####################                                        
                        if len(allele_alt_split) == 3:#A,G,<NON_REF>
                            ALT1=allele_alt_split[0]
                            ALT2=allele_alt_split[1]
                            info_field=e[8] #4 possibilite GT:AD:DP:GQ:PL, etc
                            info_field_split=info_field.split(":") #str(info_field_split[1])=="AD"
                            if str(info_field_split[0])=="GT" and str(info_field_split[1])=="AD" and str(info_field_split[2])=="DP" and str(info_field_split[3])=="GQ":
                                info_field_sample=e[9]#0/0:1:3:1:0,3,35
                                info_field_sample_split=info_field_sample.split(":")
                                GT=info_field_sample_split[0]
                                GQ=info_field_sample_split[3]
                                DP=info_field_sample_split[2]
                                if (str(GT)=="0/0" or str(GT)=="0|0") and int(GQ)>=1 and int(DP)>=1:
                                    sample_GT=str(allele_ref_new)+str(allele_ref_new)
                                    fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="0/1" or str(GT)=="0|1" or str(GT)=="1|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(allele_ref_new)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="1/1" or str(GT)=="1|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(new_ALT1)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="0/2" or str(GT)=="0|2" or str(GT)=="2|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(allele_ref_new)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="2/2" or str(GT)=="2|2") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(new_ALT2)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="1/2" or str(GT)=="1|2" or str(GT)=="2|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT2=str(ALT2)
                                            sample_GT=str(new_ALT1)+str(new_ALT2)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")

                    ####################                   
                    #####T, A, G,<NON_REF>
                    ####################

                        if len(allele_alt_split) == 4:#T,A,G,<NON_REF>
                            ALT1=allele_alt_split[0]
                            ALT2=allele_alt_split[1]
                            ALT3=allele_alt_split[2]
                            info_field=e[8] #4 possibilite GT:AD:DP:GQ:PL, etc
                            info_field_split=info_field.split(":") #str(info_field_split[1])=="AD"
                            if str(info_field_split[0])=="GT" and str(info_field_split[1])=="AD" and str(info_field_split[2])=="DP" and str(info_field_split[3])=="GQ":
                                info_field_sample=e[9]#0/0:1:3:1:0,3,35
                                info_field_sample_split=info_field_sample.split(":")
                                GT=info_field_sample_split[0]
                                GQ=info_field_sample_split[3]
                                DP=info_field_sample_split[2]
                                if (str(GT)=="0/0" or str(GT)=="0|0") and int(GQ)>=1 and int(DP)>=1:
                                    sample_GT=str(allele_ref_new)+str(allele_ref_new)
                                    fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="0/1" or str(GT)=="0|1" or str(GT)=="1|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(allele_ref_new)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="1/1" or str(GT)=="1|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(new_ALT1)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="0/2" or str(GT)=="0|2" or str(GT)=="2|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(allele_ref_new)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="2/2" or str(GT)=="2|2") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(new_ALT2)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="1/2" or str(GT)=="1|2" or str(GT)=="2|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT2=str(ALT2)
                                            sample_GT=str(new_ALT1)+str(new_ALT2)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="3/3" or str(GT)=="3|3") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                        new_ALT3=str(ALT3)
                                        sample_GT=str(new_ALT3)+str(new_ALT3)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                                 
                                if (str(GT)=="1/3" or str(GT)=="3|1" or str(GT)=="1|3") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT3=str(ALT3)
                                            sample_GT=str(new_ALT1)+str(new_ALT3)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="2/3" or str(GT)=="2|3" or str(GT)=="3|2") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                            new_ALT2=str(ALT2)
                                            new_ALT3=str(ALT3)
                                            sample_GT=str(new_ALT2)+str(new_ALT3)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="0/3" or str(GT)=="0|3" or str(GT)=="3|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                        new_ALT3=str(ALT3)
                                        sample_GT=str(allele_ref_new)+str(new_ALT3)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                              


                    ####################                   
                    ##### *, T, A, G,<NON_REF>
                    ####################

                        if len(allele_alt_split) == 5:# *,T,A,G,<NON_REF>
                            ALT1=allele_alt_split[0]
                            ALT2=allele_alt_split[1]
                            ALT3=allele_alt_split[2]
                            ALT4=allele_alt_split[3]
                            info_field=e[8] #4 possibilite GT:AD:DP:GQ:PL, etc
                            info_field_split=info_field.split(":") #str(info_field_split[1])=="AD"
                            if str(info_field_split[0])=="GT" and str(info_field_split[1])=="AD" and str(info_field_split[2])=="DP" and str(info_field_split[3])=="GQ":
                                info_field_sample=e[9]#0/0:1:3:1:0,3,35
                                info_field_sample_split=info_field_sample.split(":")
                                GT=info_field_sample_split[0]
                                GQ=info_field_sample_split[3]
                                DP=info_field_sample_split[2]
                                if (str(GT)=="0/0" or str(GT)=="0|0") and int(GQ)>=1 and int(DP)>=1:
                                    sample_GT=str(allele_ref_new)+str(allele_ref_new)
                                    fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="0/1" or str(GT)=="0|1" or str(GT)=="1|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(allele_ref_new)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                 
                                if (str(GT)=="1/1" or str(GT)=="1|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        new_ALT1=str(ALT1)
                                        sample_GT=str(new_ALT1)+str(new_ALT1)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="0/2" or str(GT)=="0|2" or str(GT)=="2|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(allele_ref_new)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="2/2" or str(GT)=="2|2") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        new_ALT2=str(ALT2)
                                        sample_GT=str(new_ALT2)+str(new_ALT2)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n") 
                                if (str(GT)=="1/2" or str(GT)=="1|2" or str(GT)=="2|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT2=str(ALT2)
                                            sample_GT=str(new_ALT1)+str(new_ALT2)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="3/3" or str(GT)=="3|3") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                        new_ALT3=str(ALT3)
                                        sample_GT=str(new_ALT3)+str(new_ALT3)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")                                 
                                if (str(GT)=="1/3" or str(GT)=="3|1" or str(GT)=="1|3") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT3=str(ALT3)
                                            sample_GT=str(new_ALT1)+str(new_ALT3)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="2/3" or str(GT)=="2|3" or str(GT)=="3|2") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                            new_ALT2=str(ALT2)
                                            new_ALT3=str(ALT3)
                                            sample_GT=str(new_ALT2)+str(new_ALT3)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                if (str(GT)=="0/3" or str(GT)=="0|3" or str(GT)=="3|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                        new_ALT3=str(ALT3)
                                        sample_GT=str(allele_ref_new)+str(new_ALT3)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
                                ### ALT4:
                                if (str(GT)=="4/4" or str(GT)=="4|4") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT4)=="A" or str(ALT4)=="T" or str(ALT4)=="G" or str(ALT4)=="C":
                                        new_ALT4=str(ALT4)
                                        sample_GT=str(new_ALT4)+str(new_ALT4)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")

                                if (str(GT)=="0/4" or str(GT)=="0|4" or str(GT)=="4|0") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT4)=="A" or str(ALT4)=="T" or str(ALT4)=="G" or str(ALT4)=="C":
                                        new_ALT4=str(ALT4)
                                        sample_GT=str(allele_ref_new)+str(new_ALT4)
                                        fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")

                                if (str(GT)=="1/4" or str(GT)=="1|4" or str(GT)=="4|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT1)=="A" or str(ALT1)=="T" or str(ALT1)=="G" or str(ALT1)=="C":
                                        if str(ALT4)=="A" or str(ALT4)=="T" or str(ALT4)=="G" or str(ALT4)=="C":
                                            new_ALT1=str(ALT1)
                                            new_ALT4=str(ALT4)
                                            sample_GT=str(new_ALT1)+str(new_ALT4)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")

                                if (str(GT)=="2/4" or str(GT)=="2|4" or str(GT)=="4|1") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT2)=="A" or str(ALT2)=="T" or str(ALT2)=="G" or str(ALT2)=="C":
                                        if str(ALT4)=="A" or str(ALT4)=="T" or str(ALT4)=="G" or str(ALT4)=="C":
                                            new_ALT2=str(ALT2)
                                            new_ALT4=str(ALT4)
                                            sample_GT=str(new_ALT2)+str(new_ALT4)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")

                                if (str(GT)=="3/4" or str(GT)=="3|4" or str(GT)=="4|3") and int(GQ)>=1 and int(DP)>=1:
                                    if str(ALT3)=="A" or str(ALT3)=="T" or str(ALT3)=="G" or str(ALT3)=="C":
                                        if str(ALT4)=="A" or str(ALT4)=="T" or str(ALT4)=="G" or str(ALT4)=="C":
                                            new_ALT3=str(ALT3)
                                            new_ALT4=str(ALT4)
                                            sample_GT=str(new_ALT3)+str(new_ALT4)
                                            fR.write(str(new_chr)+"\t"+str(new_pos)+"\t"+str(allele_ref_new)+"\t"+str(sample_GT)+"\t"+str(DP)+"\t"+str(GQ)+"\n")
