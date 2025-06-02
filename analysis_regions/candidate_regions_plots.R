#!/usr/bin/env Rscript

# Script to obtain plots for the canditate regions with the spatial dataset.
# The mean heterozygosity per sample across the candidate region was obtain with VCFtools (v. 0.1.17) : vcftools --vcf chrX_XXX_XXXMb.vcf --het --out chrx_XXX_XXMb

library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(Rmisc) #Version: 1.5.1 for summarySE() function
library(SNPRelate) #Version: 1.30.1
library(rstatix) #Version: 0.7.2

### Region chr1:

chr1_vcf="regions/spatial_dataset_chr1_100-161Mb.vcf.gz"
chr1_gds="regions/spatial_dataset_chr1_100-161Mb.gds"

snpgdsVCF2GDS(chr1_vcf, chr1_gds, method="biallelic.only")
genofile_chr1 = snpgdsOpen(chr1_gds) 
pca_result_chr1 = snpgdsPCA(genofile_chr1, autosome.only=F, remove.monosnp=T) 
tab_pca_chr1 = data.frame(sampleID = pca_result_chr1$sample.id,
                 EV1 = pca_result_chr1$eigenvect[,1],  
                 EV2 = pca_result_chr1$eigenvect[,2],   
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr1) = tab_pca_chr1$sampleID 

varPC1_chr1 = round(pca_result_chr1$varprop[1]*100,0)
varPC2_chr1 = round(pca_result_chr1$varprop[2]*100,0)

plot_PCA_chr1 <- ggplot(tab_pca_chr1, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID)) + 
  geom_point(size=4, alpha=0.5, color="blue") + xlab(paste0("PC1 (",varPC1_chr1," %)")) + ylab(paste0("PC2 (",varPC2_chr1," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") 

ggsave("img/Local_PCA_chr1.png",plot = last_plot(),width = 7, height = 5)

##############################################################################################################################

### Region chr2:

## PCA:

chr2_vcf="regions/spatial_dataset_chr2_131-186Mb.vcf.gz"
chr2_gds="regions/spatial_dataset_chr2_131-186Mb.gds"

snpgdsVCF2GDS(chr2_vcf, chr2_gds, method="biallelic.only")
genofile_chr2 = snpgdsOpen(chr2_gds) 
pca_result_chr2 = snpgdsPCA(genofile_chr2, autosome.only=F, remove.monosnp=T) 
tab_pca_chr2 = data.frame(sampleID = pca_result_chr2$sample.id,
                 EV1 = pca_result_chr2$eigenvect[,1],  
                 EV2 = pca_result_chr2$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr2) = tab_pca_chr2$sampleID 

varPC1_chr2 = round(pca_result_chr2$varprop[1]*100,0)
varPC2_chr2 = round(pca_result_chr2$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr2<-tab_pca_chr2%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr2%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr2,6)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr2,6)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr2<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr2 <- ggplot(tab_cluster_kmean_lowest_chr2, aes(x=as.numeric(EV1), y=as.numeric(-EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1_chr2," %)")) + ylab(paste0("PC2 (",varPC2_chr2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "orange", "1" = "purple", "3" = "brown", "4" = "red", "5"= "blue", "6"="black")) 

ggsave("img/Local_PCA_chr2.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr2<-read_delim("heterozygosity/spatial_dataset_chr2_131-186Mb.het") 
het_site_clusters_chr2<-left_join(het_site_chr2,tab_cluster_kmean_lowest_chr2)
het_site_clusters_chr2_het <- het_site_clusters_chr2 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr2_het_stats <- summarySE(het_site_clusters_chr2_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr2 <- ggplot(het_site_clusters_chr2_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr2_het_stats$cluster), c("6","2","3","4","1","5"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr2_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "orange", "1" = "purple", "3" = "brown", "4" = "red", "5"= "blue", "6"="black"))  

ggsave("img/Heterozygosity_chr2.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr2_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################

### Region chr3:

## PCA:

chr3_vcf="regions/spatial_dataset_chr3_101-188Mb.vcf.gz"
chr3_gds="regions/spatial_dataset_chr3_101-188Mb.gds"

snpgdsVCF2GDS(chr3_vcf, chr3_gds, method="biallelic.only")
genofile_chr3 = snpgdsOpen(chr3_gds) 
pca_result_chr3 = snpgdsPCA(genofile_chr3, autosome.only=F, remove.monosnp=T) 
tab_pca_chr3 = data.frame(sampleID = pca_result_chr3$sample.id,
                 EV1 = pca_result_chr3$eigenvect[,1],  
                 EV2 = pca_result_chr3$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr3) = tab_pca_chr3$sampleID 

varPC1 = round(pca_result_chr3$varprop[1]*100,0)
varPC2 = round(pca_result_chr3$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr3<-tab_pca_chr3%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr3%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr3,5)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr3,5)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr3<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr3 <- ggplot(tab_cluster_kmean_lowest_chr3, aes(x=as.numeric(-EV1), y=as.numeric(-EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "orange", "1" = "red", "3" = "blue", "4" = "brown", "5"= "purple", "6"="black")) 

ggsave("img/Local_PCA_chr3.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr3<-read_delim("heterozygosity/spatial_dataset_chr3_101-189Mb.het") 
het_site_clusters_chr3<-left_join(het_site_chr3,tab_cluster_kmean_lowest_chr3)
het_site_clusters_chr3_het <- het_site_clusters_chr3 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr3_het_stats <- summarySE(het_site_clusters_chr3_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr3 <- ggplot(het_site_clusters_chr3_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr3_het_stats$cluster), c("2","4","1","5","3"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr3_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "orange", "1" = "red", "3" = "blue", "4" = "brown", "5"= "purple", "6"="black"))  

ggsave("img/Heterozygosity_chr3.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr3_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################

### Region chr4:

## PCA:

chr4_vcf="regions/spatial_dataset_chr4_118-160Mb.vcf.gz"
chr4_gds="regions/spatial_dataset_chr4_118-160Mb.gds"

snpgdsVCF2GDS(chr4_vcf, chr4_gds, method="biallelic.only")
genofile_chr4 = snpgdsOpen(chr4_gds) 
pca_result_chr4 = snpgdsPCA(genofile_chr4, autosome.only=F, remove.monosnp=T) 
tab_pca_chr4 = data.frame(sampleID = pca_result_chr4$sample.id,
                 EV1 = pca_result_chr4$eigenvect[,1],  
                 EV2 = pca_result_chr4$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr4) = tab_pca_chr4$sampleID 

varPC1 = round(pca_result_chr4$varprop[1]*100,0)
varPC2 = round(pca_result_chr4$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr4<-tab_pca_chr4%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr4%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr4,2)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr4,2)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr4<-cbind(tab_pca_infos_to_keep, df)

max_PC1=max(tab_pca_kmeans_chr4$EV1)

tab_cluster_kmean_lowest_chr4 <- mutate(tab_cluster_kmean_lowest_chr4, cluster = case_when(
  EV1 == as.numeric(max_PC1) ~ 3, 
  TRUE   ~ cluster 
))

tab_cluster_kmean_lowest_chr4

plot_PCA_chr4 <- ggplot(tab_cluster_kmean_lowest_chr4, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "purple", "1" = "blue", "3" = "red")) 

ggsave("img/Local_PCA_chr4.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr4<-read_delim("heterozygosity/spatial_dataset_chr4_118-160Mb.het") 
het_site_clusters_chr4<-left_join(het_site_chr4,tab_cluster_kmean_lowest_chr4)
het_site_clusters_chr4_het <- het_site_clusters_chr4 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr4_het_stats <- summarySE(het_site_clusters_chr4_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr4 <- ggplot(het_site_clusters_chr4_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr4_het_stats$cluster), c("3", "2","1"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr4_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "purple", "1" = "blue", "3" = "red"))  

ggsave("img/Heterozygosity_chr4.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr4_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################

### Region chr6:

## PCA:

chr6_vcf="regions/spatial_dataset_chr6_132-201Mb.vcf.gz"
chr6_gds="regions/spatial_dataset_chr6_132-201Mb.gds"

snpgdsVCF2GDS(chr6_vcf, chr6_gds, method="biallelic.only")
genofile_chr6 = snpgdsOpen(chr6_gds) 
pca_result_chr6 = snpgdsPCA(genofile_chr6, autosome.only=F, remove.monosnp=T) 
tab_pca_chr6 = data.frame(sampleID = pca_result_chr6$sample.id,
                 EV1 = pca_result_chr6$eigenvect[,1],  
                 EV2 = pca_result_chr6$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr6) = tab_pca_chr6$sampleID 

varPC1 = round(pca_result_chr6$varprop[1]*100,0)
varPC2 = round(pca_result_chr6$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr6<-tab_pca_chr6%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr6%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr6,6)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr6,6)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr6<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr6 <- ggplot(tab_cluster_kmean_lowest_chr6, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "brown", "1" = "orange", "3" = "blue", "4" = "purple", "5"= "red", "6"="black")) 

ggsave("img/Local_PCA_chr6.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr6<-read_delim("heterozygosity/spatial_dataset_chr6_132-201Mb.het") 
het_site_clusters_chr6<-left_join(het_site_chr6,tab_cluster_kmean_lowest_chr6)
het_site_clusters_chr6_het <- het_site_clusters_chr6 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr6_het_stats <- summarySE(het_site_clusters_chr6_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr6 <- ggplot(het_site_clusters_chr6_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr6_het_stats$cluster), c("6","1","2","5","4","3"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr6_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "brown", "1" = "orange", "3" = "blue", "4" = "purple", "5"= "red", "6"="black"))

ggsave("img/Heterozygosity_chr6.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr6_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################

### Region chr6_2:

## PCA:

chr6_2_vcf="regions/spatial_dataset_chr6_224-229Mb.vcf.gz"
chr6_2_gds="regions/spatial_dataset_chr6_224-229Mb.gds"

snpgdsVCF2GDS(chr6_2_vcf, chr6_2_gds, method="biallelic.only")
genofile_chr6_2 = snpgdsOpen(chr6_2_gds) 
pca_result_chr6_2 = snpgdsPCA(genofile_chr6_2, autosome.only=F, remove.monosnp=T) 
tab_pca_chr6_2 = data.frame(sampleID = pca_result_chr6_2$sample.id,
                 EV1 = pca_result_chr6_2$eigenvect[,1],  
                 EV2 = pca_result_chr6_2$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr6_2) = tab_pca_chr6_2$sampleID 

varPC1 = round(pca_result_chr6_2$varprop[1]*100,0)
varPC2 = round(pca_result_chr6_2$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr6_2<-tab_pca_chr6_2%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr6_2%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr6_2,6)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr6_2,6)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr6_2<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr6_2 <- ggplot(tab_cluster_kmean_lowest_chr6_2, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "blue", "1" = "brown", "3" = "orange", "4" = "purple", "5"= "red", "6"="black")) 

ggsave("img/Local_PCA_chr6_2.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr6_2<-read_delim("heterozygosity/spatial_dataset_chr6_224-229Mb.het") 
het_site_clusters_chr6_2<-left_join(het_site_chr6_2,tab_cluster_kmean_lowest_chr6_2)
het_site_clusters_chr6_2_het <- het_site_clusters_chr6_2 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr6_2_het_stats <- summarySE(het_site_clusters_chr6_2_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr6_2 <- ggplot(het_site_clusters_chr6_2_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr6_2_het_stats$cluster), c("6","3","1","5","4","2"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr6_2_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "blue", "1" = "brown", "3" = "orange", "4" = "purple", "5"= "red", "6"="black")) 

ggsave("img/Heterozygosity_chr6_2.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr6_2_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################

### Region chr7:

## PCA:

chr7_vcf="regions/spatial_dataset_chr7_128-177Mb.vcf.gz"
chr7_gds="regions/spatial_dataset_chr7_128-177Mb.gds"

snpgdsVCF2GDS(chr7_vcf, chr7_gds, method="biallelic.only")
genofile_chr7 = snpgdsOpen(chr7_gds) 
pca_result_chr7 = snpgdsPCA(genofile_chr7, autosome.only=F, remove.monosnp=T) 
tab_pca_chr7 = data.frame(sampleID = pca_result_chr7$sample.id,
                 EV1 = pca_result_chr7$eigenvect[,1],  
                 EV2 = pca_result_chr7$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr7) = tab_pca_chr7$sampleID 

varPC1 = round(pca_result_chr7$varprop[1]*100,0)
varPC2 = round(pca_result_chr7$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr7<-tab_pca_chr7%>%select(EV1)
tab_pca_infos_to_keep<-tab_pca_chr7%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr7,3)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr7,3)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr7<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr7 <- ggplot(tab_cluster_kmean_lowest_chr7, aes(x=as.numeric(-EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "blue", "1" = "purple", "3" = "red")) 

ggsave("img/Local_PCA_chr7.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr7<-read_delim("heterozygosity/spatial_dataset_chr7_128-177Mb.het") 
het_site_clusters_chr7<-left_join(het_site_chr7,tab_cluster_kmean_lowest_chr7)
het_site_clusters_chr7_het <- het_site_clusters_chr7 %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr7_het_stats <- summarySE(het_site_clusters_chr7_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr7 <- ggplot(het_site_clusters_chr7_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr7_het_stats$cluster), c("3","1","2"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr7_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") + scale_color_manual(values = c("2"= "blue", "1" = "purple", "3" = "red"))   

ggsave("img/Heterozygosity_chr7.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr7_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test

##############################################################################################################################
##############################################################################################################################

### Region chr3 temporal dataset:

## PCA:

chr3_temporal_vcf="regions/temporal_dataset_chr3_101-189Mb.vcf.gz"
chr3_temporal_gds="regions/temporal_dataset_chr3_101-189Mb.gds"

snpgdsVCF2GDS(chr3_temporal_vcf, chr3_temporal_gds, method="biallelic.only")
genofile_chr3_temporal = snpgdsOpen(chr3_temporal_gds) 
pca_result_chr3_temporal = snpgdsPCA(genofile_chr3_temporal, autosome.only=F, remove.monosnp=T) 
tab_pca_chr3_temporal = data.frame(sampleID = pca_result_chr3_temporal$sample.id,
                 EV1 = pca_result_chr3_temporal$eigenvect[,1],  
                 EV2 = pca_result_chr3_temporal$eigenvect[,2],  
                 stringsASFactors=FALSE 
)

rownames(tab_pca_chr3_temporal) = tab_pca_chr3_temporal$sampleID 

varPC1 = round(pca_result_chr3_temporal$varprop[1]*100,0)
varPC2 = round(pca_result_chr3_temporal$varprop[2]*100,0)

## Kmean clustering and plot of PCA:

tab_pca_kmeans_chr3_temporal<-tab_pca_chr3_temporal%>%select(EV1,EV2)
tab_pca_infos_to_keep<-tab_pca_chr3_temporal%>%select(sampleID,EV1,EV2)

set.seed(123)
lowest_r<-kmeans(tab_pca_kmeans_chr3_temporal,6)
lowest_totwithinss<-lowest_r$tot.withinss
for (i in 1:1000) {
  r_kmean<-kmeans(tab_pca_kmeans_chr3_temporal,6)
  new_totwithinss<-r_kmean$tot.withinss
  if (as.numeric(new_totwithinss)<as.numeric(lowest_totwithinss)) {
    lowest_totwithinss<-new_totwithinss
    lowest_r<-r_kmean
  }
}

df<-data.frame(lowest_r["cluster"])
tab_cluster_kmean_lowest_chr3_temporal<-cbind(tab_pca_infos_to_keep, df)

plot_PCA_chr3_temporal <- ggplot(tab_cluster_kmean_lowest_chr3_temporal, aes(x=as.numeric(EV1), y=as.numeric(EV2), label=sampleID, col=as.character(cluster))) + 
  geom_point(size=4, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.title = element_text(size=18),
        axis.title.x = element_text(size=18,vjust=-0.4),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18) 
    ) +labs(col = "") + theme(legend.position="none") +  scale_color_manual(values = c("2"= "orange", "1" = "blue", "3" = "brown", "4" = "purple", "5"= "red", "6"="black"))  

ggsave("img/Local_PCA_chr3_temporal.png",plot = last_plot(),width = 7, height = 5)

## Heterozygosity

het_site_chr3_temporal<-read_delim("heterozygosity/temporal_dataset_chr3_101-189Mb.het") 
het_site_clusters_chr3_temporal<-left_join(het_site_chr3_temporal,tab_cluster_kmean_lowest_chr3_temporal)
het_site_clusters_chr3_temporal_het <- het_site_clusters_chr3_temporal %>%
  select(sampleID, HOM, N_SITES, cluster, EV1, EV2) %>%
  mutate(HET = (N_SITES-HOM)/N_SITES, na.rm = TRUE) 
het_site_clusters_chr3_temporal_het_stats <- summarySE(het_site_clusters_chr3_temporal_het, measurevar="HET", groupvars=c("cluster"))

het_plot_chr3_temporal <- ggplot(het_site_clusters_chr3_temporal_het_stats, aes(x=fct_rev(fct_relevel(as.character(het_site_clusters_chr3_temporal_het_stats$cluster), c("6","2","3","5","4","1"))), y=HET,col=as.character(cluster)))  +
    geom_errorbar(width=.2, aes(ymin=HET-(1.96*se), ymax=HET+1.96*se)) +
    geom_point(size=5, aes(color = as.character(het_site_clusters_chr3_temporal_het_stats$cluster))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(size=18),
        strip.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18,vjust=1.5)) + xlab(" ") + ylab("Heterozygosity") + theme(legend.position = "none") +  scale_color_manual(values = c("2"= "orange", "1" = "blue", "3" = "brown", "4" = "purple", "5"= "red", "6"="black"))  

ggsave("img/Heterozygosity_chr3_temporal_dataset.png",plot = last_plot(),width = 7, height = 5)

stat.test <- het_site_clusters_chr3_temporal_het %>% 
  wilcox_test(HET ~ cluster) %>%
  add_significance()
stat.test





