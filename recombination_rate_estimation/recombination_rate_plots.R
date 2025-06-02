#!/usr/bin/env Rscript

# Script to obtain plots of the estimated recombination rate across the chromosomes. 
# Estimation of te recombination rate per non-overlapping windows across the chromosomes was obtained with ReLERNN (v. 1.0.O): https://github.com/kr-colab/ReLERNN 
# We included the windows in the candidate regions when more than half of their length covered a candidate regins previously identified with a population genomic approach.
# For this reason, the candidate region defined here can begin and finish a few Mbases before or after the positions identified previously with a population genomic approach. 
# Can be lauched as: ./recombination_rate_plots.R -c chr2 -s 126016000 -e 189024000

library(optparse) #Version: 1.7.3
library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(scales) #Version: 1.2.1

option_list = list(
  make_option(c("-c","--chr"), action="store", type="character", help="chromosome"),
  make_option(c("-s","--start"), action="store", type="integer", help="start position of the region to plot"),
  make_option(c("-e","--end"), action="store", type="integer", help="end position of the region to plot")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

chromosome = opt$chr
startP = opt$start
endP = opt$end

file_all<-"relernn_output/spatial_dataset_relernn.PREDICT.BSCORRECTED.txt"

recomb_rate_all <- read.delim(file_all, sep="\t")

tab_recomb<-recomb_rate_all%>%
  mutate(variable = case_when(
      (chrom == "b'chr1'" & start >= 101586000 & end <= 169310000) ~ "Chr1 99-162Mb",
      (chrom == "b'chr2'" & start >= 126016000 & end <= 189024000) ~ "Chr2 131-186Mb",
      (chrom == "b'chr3'" & start >= 101238580 & end <= 189497000) ~ "Chr3 101-189Mb",
      (chrom == "b'chr4'" & start >= 116754000 & end <= 160157576) ~ "Chr4 118-160Mb",
      (chrom == "b'chr6'" & start >= 132717319 & end <= 209256000) ~ "Chr6 132-201Mb",
      (chrom == "b'chr7'" & start >= 128666580 & end <= 178080000) ~ "Chr7 128-177Mb"))

tab_recomb_var <- tab_recomb$variable %>% replace_na("All genome")
tab_recomb$zone <- tab_recomb_var

region_chr1<-tab_recomb%>%filter(zone=="Chr1 99-162Mb") 
region_chr2<-tab_recomb%>%filter(zone=="Chr2 131-186Mb")
region_chr3<-tab_recomb%>%filter(zone=="Chr3 101-189Mb")
region_chr4<-tab_recomb%>%filter(zone=="Chr4 118-160Mb")
region_chr6<-tab_recomb%>%filter(zone=="Chr6 132-201Mb")
region_chr7<-tab_recomb%>%filter(zone=="Chr7 128-177Mb")
region_all<-tab_recomb%>%filter(zone=="All genome")

wilcox.test(region_chr1$recombRate,recomb_rate_all$recombRate)
wilcox.test(region_chr2$recombRate,recomb_rate_all$recombRate)
wilcox.test(region_chr3$recombRate,recomb_rate_all$recombRate)
wilcox.test(region_chr4$recombRate,recomb_rate_all$recombRate)
wilcox.test(region_chr6$recombRate,recomb_rate_all$recombRate)
wilcox.test(region_chr7$recombRate,recomb_rate_all$recombRate)

chr<-recomb_rate_all%>%filter(chrom==paste0("b'",chromosome,"'"))

plot<-ggplot(chr, aes(x=start, y=recombRate, ymax = CI95HI, ymin = CI95LO)) + annotate('rect', xmin=as.numeric(startP), xmax=as.numeric(endP), ymin=0, ymax=0.000000011, alpha=.2, fill='orange') +  geom_line(size=1, alpha=0.5) + xlab(paste0(chromosome)) + ylab("Recombination rate estimation") + geom_ribbon(alpha=0.1) +
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.y = element_text(size=15),
        axis.text.x= element_text(size=20,angle=75,hjust=1,vjust=1),
        strip.text.x = element_text(size = 25),
        axis.title.x = element_text(size=23,vjust=-0.4),
        axis.title.y = element_text(size=23,vjust=1.5)
    ) + scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6), breaks = c(50000000,100000000,150000000,200000000,250000000,300000000)) 


ggsave(paste0("img/",chromosome,"_estimation_recombination_rate.png"), width = 50, height = 20, units = "cm")


















