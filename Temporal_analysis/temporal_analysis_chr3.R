#!/usr/bin/env Rscript

# Script to obtain plots for the temporal analysis of the candidate regions of chromosome 3 with the temporal dataset.
# The mean heterozygosity per sample across the candidate region was obtain with VCFtools (v. 0.1.17) : vcftools --vcf chr3_101_189Mb_temporal.vcf --het --out chr3_101_189Mb_temporal

library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(Rmisc) #Version: 1.5.1 for summarySE() function
library(SNPRelate) #Version: 1.30.1
library(rstatix) #Version: 0.7.2

## PCA:

chr3_vcf="temporal_dataset_chr3_101-189Mb.vcf.gz"
chr3_gds="temporal_dataset_chr3_101-189Mb.gds"

snpgdsVCF2GDS(chr3_vcf, chr3_gds, method="biallelic.only")
genofile_chr3 = snpgdsOpen(chr3_gds) 
pca_result_chr3 = snpgdsPCA(genofile_chr3, autosome.only=F, remove.monosnp=T) 
tab_pca_chr3 = data.frame(sampleID = pca_result_chr3$sample.id,
                 EV1 = pca_result_chr3$eigenvect[,1],  
                 EV2 = pca_result_chr3$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr3) = tab_pca_chr3$sampleID 

varPC1_chr3 = round(pca_result_chr3$varprop[1]*100,0)
varPC2_chr3 = round(pca_result_chr3$varprop[2]*100,0)

## PCA with the years:

infos_temp<- read.delim("temporal_data.txt")
tab_pca_infos_temp<-left_join(tab_pca_chr3,infos_temp)

acp_temp<-ggplot(tab_pca_infos_temp, aes(x=EV1, y=EV2, label=sampleID, col=as.character(year), shape=pop)) + 
  geom_point(size=3, alpha=0.5) + xlab(paste0("PC1 (",varPC1_chr3," %)")) + ylab(paste0("PC2 (",varPC2_chr3," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=14), 
        axis.text = element_text(size=14),
        legend.text = element_text(size=14)
    ) +labs(col = "") 
acp_temp

ggsave("img/Local_PCA_chr3_temporal.png",width = 7, height = 5)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr3<-tab_pca_chr3%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_infos_temp%>%select(sampleID,EV1,EV2,year)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr3,6)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr3,6)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr3<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr3 <- ggplot(tab_cluster_kmean_lowest_chr3, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1_chr3," %)")) + ylab(paste0("PC2 (",varPC2_chr3," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "orange", "1" = "blue", "3" = "brown", "4" = "purple", "5"= "red", "6"="black")) 

ggsave("img/Local_PCA_chr3_temporal_kmean_clusters.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr3<-read_delim("chr3_101_189Mb_temporal.het") 
het_site_clusters_chr3<-left_join(het_site_chr3,tab_cluster_kmean_lowest_chr3)
het_site_clusters_chr3_het <- het_site_clusters_chr3 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr3_het_stats <- summarySE(het_site_clusters_chr3_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr3 <- ggplot(het_site_clusters_chr3_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr3_het_stats$cluster), c("6","2","3","5","4","1"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr3_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") +  scale_color_manual(values = c("2"= "orange", "1" = "blue", "3" = "brown", "4" = "purple", "5"= "red", "6"="black"))  

ggsave("img/Heterozygosity_chr3.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr3_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

## Fischer tests:

print("Fischer tests:")

# 2016
year2016<-tab_cluster_kmean_lowest_chr3%>%filter(year=="2016")

cluster1_2016<-year2016%>%filter(cluster=="1") #blue AA
AA_2016<-as.numeric(length(cluster1_2016$sampleID))

cluster2_2016<-year2016%>%filter(cluster=="4") #purple AB
AB_2016<-as.numeric(length(cluster2_2016$sampleID))

cluster3_2016<-year2016%>%filter(cluster=="5") #red BB
BB_2016<-as.numeric(length(cluster3_2016$sampleID))

cluster4_2016<-year2016%>%filter(cluster=="3") #brown AC
AC_2016<-as.numeric(length(cluster4_2016$sampleID))

cluster5_2016<-year2016%>%filter(cluster=="2") #orange BC
BC_2016<-as.numeric(length(cluster5_2016$sampleID))

cluster6_2016<-year2016%>%filter(cluster=="6") #black CC
CC_2016<-as.numeric(length(cluster6_2016$sampleID))

# 1976
year1976<-tab_cluster_kmean_lowest_chr3%>%filter(year=="1976")

cluster1_1976<-year1976%>%filter(cluster=="1") #blue AA
AA_1976<-as.numeric(length(cluster1_1976$sampleID))

cluster2_1976<-year1976%>%filter(cluster=="4") #purple AB
AB_1976<-as.numeric(length(cluster2_1976$sampleID))

cluster3_1976<-year1976%>%filter(cluster=="5") #red BB
BB_1976<-as.numeric(length(cluster3_1976$sampleID))

cluster4_1976<-year1976%>%filter(cluster=="3") #brown AC
AC_1976<-as.numeric(length(cluster4_1976$sampleID))

cluster5_1976<-year1976%>%filter(cluster=="2") #orange BC
BC_1976<-as.numeric(length(cluster5_1976$sampleID))

cluster6_1976<-year1976%>%filter(cluster=="6") #black CC
CC_1976<-as.numeric(length(cluster6_1976$sampleID))

# Overall Fischer test:

fisher.test(matrix(c(as.numeric(AA_2016),as.numeric(AA_1976),as.numeric(AB_2016),as.numeric(AB_1976),as.numeric(BB_2016),as.numeric(BB_1976),as.numeric(AC_2016),as.numeric(AC_1976),as.numeric(BC_2016),as.numeric(BC_1976),as.numeric(CC_2016),as.numeric(CC_1976)),6,2, byrow=TRUE))

display_matrix<-matrix(c(as.numeric(AA_2016),as.numeric(AA_1976),as.numeric(AB_2016),as.numeric(AB_1976),as.numeric(BB_2016),as.numeric(BB_1976),as.numeric(AC_2016),as.numeric(AC_1976),as.numeric(BC_2016),as.numeric(BC_1976),as.numeric(CC_2016),as.numeric(CC_1976)),6,2, byrow=TRUE)

print(display_matrix)

#     [,1] [,2]
#[1,]   59   61
#[2,]   59   69
#[3,]   11   23
#[4,]   40   26
#[5,]   17    9
#[6,]    3    0


samples_2016<-nrow(year2016)
samples_1976<-nrow(year1976)

## AA
no_AA_2016<-as.numeric(samples_2016)-as.numeric(AA_2016)
no_AA_1976<-as.numeric(samples_1976)-as.numeric(AA_1976)
matrix_AA<-matrix(c(as.numeric(AA_2016),as.numeric(AA_1976),as.numeric(no_AA_2016),as.numeric(no_AA_1976)), 2,2, byrow=TRUE)
print(matrix_AA)
#     [,1] [,2]
#[1,]   59   61
#[2,]  130  127

fisher.test(matrix_AA)

## AB
no_AB_2016<-as.numeric(samples_2016)-as.numeric(AB_2016)
no_AB_1976<-as.numeric(samples_1976)-as.numeric(AB_1976)
matrix_AB<-matrix(c(as.numeric(AB_2016),as.numeric(AB_1976),as.numeric(no_AB_2016),as.numeric(no_AB_1976)), 2,2, byrow=TRUE)
print(matrix_AB)
fisher.test(matrix_AB)

## BB
no_BB_2016<-as.numeric(samples_2016)-as.numeric(BB_2016)
no_BB_1976<-as.numeric(samples_1976)-as.numeric(BB_1976)
matrix_BB<-matrix(c(as.numeric(BB_2016),as.numeric(BB_1976),as.numeric(no_BB_2016),as.numeric(no_BB_1976)), 2,2, byrow=TRUE)
print(matrix_BB)
fisher.test(matrix_BB)
 

## BC
no_BC_2016<-as.numeric(samples_2016)-as.numeric(BC_2016)
no_BC_1976<-as.numeric(samples_1976)-as.numeric(BC_1976)
matrix_BC<-matrix(c(as.numeric(BC_2016),as.numeric(BC_1976),as.numeric(no_BC_2016),as.numeric(no_BC_1976)), 2,2, byrow=TRUE)
print(matrix_BC)
fisher.test(matrix_BC)
 
 







