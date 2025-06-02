#!/usr/bin/env Rscript

library(lostruct) #Version: 0.0.0.9000
library(FDRestimation) #Version: 1.0.1
library(DescTools) #Version: 0.99.49
library(tidyverse) #Version: 2.0.0
library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 1.1.2
library(scales) #Version: 1.2.1

# Script to detect genomic windows with a deviant population structure using the R package lostruct (Li & Ralph 2019, https://github.com/petrelharp/local_pca)
# The sript also allows to cluster the genomic outlier windows and to define candidate regions for putative inversions detection 

# The VCF file has first to be converted in bcf format and indexed with 
# bcftools convert -O b mydata.vcf > mydata.bcf
# bcftools index mydata.bcf

vcf_file <- "spatial_dataset.bcf"
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

# Function for clustering of the outliers genomic windows:

clustering<-function(list_parameter, gap, lenght_list) {
    c<-0
    size_list_parameter<-length(list_parameter)
    list_result <- list()
    new_l<-list()
    previous_ID<-1 
    for (i in list_parameter) { 
        c<-c+1
        ID<-i
        diff<-ID-previous_ID
        if ((diff<gap)) { 
            if (!(previous_ID %in% new_l)) {
            new_l[[(length(new_l) + 1)]] <- previous_ID
            }
            new_l[[(length(new_l) + 1)]] <- ID
        } else {
          if (length(new_l)>=lenght_list) { 
            list_result[[(length(list_result) + 1)]] <- new_l 
          }
        new_l<-list()
        }
      previous_ID<-ID
    }
    if ((length(new_l)>=lenght_list)&(diff<gap)&(c==size_list_parameter)) { 
        list_result[[(length(list_result) + 1)]] <- new_l 
    }
    return(list_result) 
}

# Identification of the candidate regions. 
# This step aims at clustering the outlier genomic windows into candidate regions. 
# Briefly, for each dimension:
# 1. We identify the chromosomes with a minimum number of outlier windows (parameter min_outlier_windows) and which gather a greater number of outlier windows compared to the rest of the genome using a G-test.
# 2. We cluster the outlier windows by allowing a configurable number (parameter nb_windows_gap) of non-outlier windows between two consecutive outliers windows.
# 3. A candidate region is defined by the start position and the end position of the first and the last genomic windows of a cluster. 
# We display candidate regions which harbor a minimum number of outlier windows (parameter nb_windows_min_cluster).

# Configurable parameters:
min_outlier_windows<- 3
Gtest_threshold<-0.05
nb_windows_gap <- 3
nb_windows_min_cluster<-3

tab_regions <- data.frame(ID=character(),MDS=character(),Chr=character(),Start=integer(),End=integer(),NbOutliers=integer(),stringsAsFactors=FALSE)
windows_tot<-nrow(tab) 
listOutliersRegions=list()
id_regions<-1 

for (i in 1:dimensions) {
  MDS<-paste0("MDS",i)
  FDR_MDS<-paste0("reject_fdr_MDS",i)
  for (j in 1:7) { #for each chromosome
    chromosome<-paste0("chr",j)
    windows_chr<-nrow(tab%>%filter(Chr==chromosome)) 
    nb_outlier_windows_chr<-nrow(tab %>% filter(tab[, FDR_MDS] == "Reject.H0") %>% filter(Chr==chromosome))
    if (nb_outlier_windows_chr > min_outlier_windows) {
      nb_outlier_w_tot<-nrow(tab %>% filter(tab[, FDR_MDS] == "Reject.H0"))
      matrix = (paste0("
           windows  significative unsignificative
           tot  ",nb_outlier_w_tot," ",windows_tot-nb_outlier_w_tot,"
           chr  ",nb_outlier_windows_chr," ",windows_chr-nb_outlier_windows_chr,"
           "))
      matrix_gtest = as.matrix(read.table(textConnection(matrix),
                         header=TRUE,
                         row.names=1))
      r_GTest<-GTest(matrix_gtest)
      G<-r_GTest[1]
      pvalue_Gtest<-r_GTest[3]
      if (pvalue_Gtest<Gtest_threshold) { 
        tab_1_p<-tab%>%select(ID,FDR_MDS,MDS,Chr,posB,posE)
        tab_2_p<-tab_1_p%>%filter(Chr==chromosome)
        tab_3_p<-tab_2_p%>%filter(tab_2_p[, MDS]<0)
        tab_4_p<-tab_3_p%>%filter(tab_3_p[, FDR_MDS] == "Reject.H0")
        if (length(tab_4_p[,1])>0) {
            list_cluster_pos<-clustering(tab_4_p[,1], nb_windows_gap, nb_windows_min_cluster) #clustering function written above
            if (length(list_cluster_pos) > 0) { 
              for (cluster_p in list_cluster_pos) { 
                nb_of_positive_outlier_windows<-length(cluster_p) 
                firstID<-cluster_p[1]
                lastID<-cluster_p[nb_of_positive_outlier_windows]
                first_w<-tab_4_p%>%filter(ID==firstID)%>%select(posB)
                first_pos<-first_w[1,1] 
                last_w<-tab_4_p%>%filter(ID==lastID)%>%select(posE)
                last_pos<-last_w[1,1]
                ID_region<-paste0(id_regions)
                tab_regions[nrow(tab_regions) + 1,] = c(ID_region,MDS,chromosome,first_pos,last_pos,nb_of_positive_outlier_windows)
                id_regions=id_regions+1
                listOutliersRegions[[(length(listOutliersRegions) + 1)]] <- cluster_p
              }
            }
        }
        tab_1_n<-tab%>%select(ID,FDR_MDS,MDS,Chr,posB,posE)
        tab_2_n<-tab_1_n%>%filter(Chr==chromosome)
        tab_3_n<-tab_2_n%>%filter(tab_2_n[, MDS]>0)
        tab_4_n<-tab_3_n%>%filter(tab_3_n[, FDR_MDS] == "Reject.H0")
        if (length(tab_4_n[,1])>0) {
            list_cluster_neg<-clustering(tab_4_n[,1], nb_windows_gap, nb_windows_min_cluster) #clustering function written above
            if (length(list_cluster_neg) > 0) { 
              for (cluster_n in list_cluster_neg) { 
                nb_of_negative_outlier_windows<-length(cluster_n)
                firstID<-cluster_n[1]
                lastID<-cluster_n[nb_of_negative_outlier_windows]
                first_w<-tab_4_n%>%filter(ID==firstID)%>%select(posB)
                first_pos<-first_w[1,1] 
                last_w<-tab_4_n%>%filter(ID==lastID)%>%select(posE)
                last_pos<-last_w[1,1] 
                ID_region<-paste0(id_regions) 
                tab_regions[nrow(tab_regions) + 1,] = c(ID_region,MDS,chromosome,first_pos,last_pos,nb_of_negative_outlier_windows)
                id_regions=id_regions+1
                listOutliersRegions[[(length(listOutliersRegions) + 1)]] <- cluster_n
              }
            }
        }
        }
      }
  }
}
tab_regions_size<-tab_regions%>%
  mutate(SizeRegions=sprintf("%.1f",(as.numeric(End)-as.numeric(Start))/1000000))
tab_regions_size_chr<-tab_regions_size%>%filter(Chr=="chr3")
tab_regions_size

write.csv(tab_regions_size, "candidate_regions_spatial_dataset.csv", row.names=TRUE)
























