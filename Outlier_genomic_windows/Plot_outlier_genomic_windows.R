#!/usr/bin/env Rscript

library(lostruct) #Version: 0.0.0.9000
library(FDRestimation) #Version: 1.0.1
library(DescTools) #Version: 0.99.49
library(tidyverse) #Version: 2.0.0
library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 1.1.2
library(scales) #Version: 1.2.1

# Script to detect and plot genomic windows with a deviant population structure using the R package lostruct (Li & Ralph 2019, https://github.com/petrelharp/local_pca). 

# The VCF file has first to be converted in bcf format and indexed with 
# bcftools convert -O b mydata.vcf > mydata.bcf
# bcftools index mydata.bcf

vcf_file <- "spatial_dataset_126_samples.bcf"
windows_size <- 100
principal_components <- 2 
dimensions <- 40 
FDR_threshold <- 0.05

# Briefly, the functions of the R package lostruct: 
#1) divides the genome in non ovelapping windows composed of 100 SNPs, 
#2) performs a local PCA with each window, 
#3) computes a pairwise comparison matrix of the local PCAs, 
#4) performs a multidimensionnal scaling analysis based on the matrix of comparison to differentiate the genomic windows.

snps <- vcf_windower(vcf_file,size=windows_size,type='snp')
pcs <- eigen_windows(snps,k=principal_components)
pcdist <- pc_dist(pcs,npc=principal_components)
nas <- is.na(pcdist[,1])
fit40d<-cmdscale(pcdist[!nas,!nas],eig=TRUE,k=40)

# The python script get_windows_positions.py has to be used to obtain positions of the windows, as follows:
# get_windows_positions.py -v VCF_FILE (spatial_dataset.vcf) -w WINDOWS_SIZE (100 )-o OUTPUT (windows_positions_spatial_dataset.txt)
# and allows to obtain the file: windows_positions_spatial_dataset.txt

tab <- read.delim("windows_positions_spatial_dataset.txt")

for (i in 1:dimensions) {
  MDS<-fit40d$points[,i]
  df_coord <- as.data.frame(MDS) 
  tab<-cbind(df_coord,tab)
  names(tab)[names(tab) == 'MDS'] <- c(paste0("MDS", i))
}

# The values of the multidimensional scaling analysis are centered and reduced for each dimension:
for (i in 1:dimensions) {
  tab<-tab %>% 
    mutate(MDS_c_r=(tab[, paste0("MDS", i)]-mean(tab[, paste0("MDS", i)]))/sd(tab[, paste0("MDS", i)]))
  names(tab)[names(tab) == 'MDS_c_r'] <- c(paste0("MDS",i,"_c_r"))
}
for (i in 1:dimensions) {
  tab<-tab %>% 
    mutate(p_value_MDS=pnorm(abs(tab[, paste0("MDS",i,"_c_r")]), mean = 0, sd = 1, lower.tail = FALSE))
  names(tab)[names(tab) == 'p_value_MDS'] <- c(paste0("p_value_MDS",i))
}

# Identification of outliers values for each dimension of the multidimensional scaling analysis with a false discovery rate correction of 5%:
for (i in 1:dimensions) { 
  fdr_res_MDS<-p.fdr(pvalues=tab[, paste0("p_value_MDS",i)],threshold=FDR_threshold)
  reject_fdr_MDS<-fdr_res_MDS$`Reject Vector`
  df_reject_fdr_MDS <- as.data.frame(reject_fdr_MDS) 
  tab <- cbind(tab, df_reject_fdr_MDS)
  names(tab)[names(tab) == 'reject_fdr_MDS'] <- c(paste0("reject_fdr_MDS", i))
}

# Plot of the values of the multidimensional scaling analysis:
for (i in 1:dimensions) { 
  MDS<-paste0("MDS",i)
  print(MDS)
  reject<-paste0("reject_fdr_MDS",i)
  plot_outliers<-ggplot(tab, aes(x=posB, y=tab[, paste0("MDS",i)], col=tab[, paste0("reject_fdr_MDS",i)])) +  
    geom_point(size=2, alpha=0.5) + xlab("") + ylab(paste0("MDS",i)) +
    theme(legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
          axis.text = element_text(size=12),
          axis.text.x= element_text(size=8,angle=75,hjust=1,vjust=1),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_text(size=12,vjust=-0.4),
          axis.title.y = element_text(size=12,vjust=1.5)
      ) +       
    facet_grid (.~ Chr, scales = "free_x", space = "free_x") +
    scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6), breaks = c(50000000,100000000,150000000,200000000,250000000,300000000)) + theme(legend.title=element_blank())  +  scale_colour_discrete(labels=c("no outliers","outliers")) + scale_color_manual(values = c("Reject.H0" = "red","FTR.H0" = "blue")) 
  ggsave(paste0("img_outliers_genomic_windows/MDS",i,".png"),width=7, height=4)
}


