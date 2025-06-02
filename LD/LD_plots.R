#!/usr/bin/env Rscript

# Script to obtain plot of r² values across the chromosomes.
# We first reduced the total number of SNPs using VCFtools (v. 0.1.16) and the--thin parameter: vcftools --vcf spatial_dataset_chrX.vcf --thin 1000 --recode --recode-INFO-all --out spatial_dataset_chrX_thin1000.vcf
# The r² values were then obtained with plink (v. 1.9) for each pair of SNPs: plink --vcf spatial_dataset_chrX_thin1000.vcf --r2 –ld-window-r2 0.0 --ld-window-kb 310000000 --ld-window 1000000 --out spatial_dataset_chrX_thin1000_plink.ld

library(optparse) #Version: 1.7.3
library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(scales) #Version: 1.2.1

option_list = list(
  make_option(c("-r","--rplink"), action="store", type="character", help="Input file with r² values"),
  make_option(c("-c","--chr"), action="store", type="character", help="chromosome")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

r_file = opt$rplink
chromosome = opt$chr

file_ld <- read.delim(r_file, sep=" ")

X=file_ld$BP_A
Y=file_ld$BP_B
value=file_ld$R2
df_plink<-data.frame(X,Y,value)

ld_plot <- ggplot(data = df_plink, aes(x=Y, y=X, fill=value)) + 
  geom_tile(width=4000000,height=4000000) +
  scale_fill_gradient(low = "white", high = "black") +
    xlab(paste0(chromosome," (Mb)")) + 
    ylab(paste0(chromosome," (Mb)")) +
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), breaks = c(50000000,100000000,150000000,200000000, 250000000)) +
  scale_y_continuous(labels = unit_format(unit = "", scale = 1e-6), breaks = c(50000000,100000000,150000000,200000000, 250000000)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.y= element_text(size=14),
        axis.text.x= element_text(size=14),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title.y = element_text(size=20,vjust=1.5),
        axis.title.x = element_text(size=20,vjust=0.5)) + labs(fill = "r²")

ggsave(paste0("img/",chromosome,"_r2.png"),plot = last_plot(),width = 7, height = 5)

