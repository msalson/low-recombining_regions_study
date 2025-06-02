#!/usr/bin/env Rscript

# Script to obtain plot from Fis values across the genome.

library(optparse) #Version: 1.7.3
library(ggplot2) #Version: 3.4.3
library(dplyr) #Version: 2.3.3
library(tidyverse) #Version: 2.0.0
library(Rmisc) #Version: 1.5.1

option_list = list(
  make_option(c("-i","--input"), action="store", type="character", help="Input file with Fis values")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fis_file = opt$input

fis<-read.delim(fis_file)

#

tab_f<-fis%>%
  mutate(variable = case_when(
      (chromosome == "chr1" & locus >= 99793192 & locus <= 161758608) ~ "Chr1 100-162 Mb",
      (chromosome == "chr2" & locus >= 131088113 & locus <= 186649394) ~ "Chr2 131-187 Mb",
      (chromosome == "chr3" & locus >= 101238580 & locus <= 188834466) ~ "Chr3 101-189 Mb",
      (chromosome == "chr4" & locus >= 118491302 & locus <= 160157576) ~ "Chr4 118-160 Mb",
      (chromosome == "chr6" & locus >= 132717319 & locus <= 201482742) ~ "Chr6 133-201 Mb",
      (chromosome == "chr6" & locus >= 224133144 & locus <= 229409049) ~ "Chr6 224-229Mb",
      (chromosome == "chr7" & locus >= 128666580 & locus <= 177591959) ~ "Chr7 129-178Mb"))

tab_variable <- tab_f$variable %>% replace_na("All genome")

tab_f$zone <- tab_variable

tab_f_stats <- summarySE(data = tab_f, measurevar="Fis", groupvars=c("zone"), na.rm = FALSE, conf.interval = 0.95, .drop = TRUE)

#####################

ggplot(tab_f_stats, aes(x=zone, y=Fis))  +
    geom_errorbar(width=.1, aes(ymin=Fis-(1.96*se), ymax=Fis+1.96*se)) +
    geom_point(shape=21, size=1, fill="black") +  geom_hline(yintercept=0, linetype="dashed", color = "gray", size=0.4)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
        axis.text.x= element_text(size=8,angle=25,hjust=1,vjust=1),
        axis.text.y= element_text(size=8),
        strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=10,vjust=1.5)) + ylab(expression(paste(F[is],sep=""))) +
    scale_y_continuous(limits=c(-0.15, 0.15), breaks=c(-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15))

#####################

tab_all <- tab_f%>%filter(zone=="All genome")
tab_chr1 <- tab_f%>%filter(zone=="Chr1 100-162 Mb")
tab_chr2 <- tab_f%>%filter(zone=="Chr2 131-187 Mb")
tab_chr3 <- tab_f%>%filter(zone=="Chr3 101-189 Mb")
tab_chr4 <- tab_f%>%filter(zone=="Chr4 118-160 Mb")
tab_chr6_1 <- tab_f%>%filter(zone=="Chr6 133-201 Mb")
tab_chr6_2 <- tab_f%>%filter(zone=="Chr6 224-229Mb")
tab_chr7 <- tab_f%>%filter(zone=="Chr7 129-178Mb")

wilcox.test(tab_all$Fis, tab_chr1$Fis)
wilcox.test(tab_all$Fis, tab_chr2$Fis)
wilcox.test(tab_all$Fis, tab_chr3$Fis)
wilcox.test(tab_all$Fis, tab_chr4$Fis)
wilcox.test(tab_all$Fis, tab_chr6_1$Fis)
wilcox.test(tab_all$Fis, tab_chr6_2$Fis)
wilcox.test(tab_all$Fis, tab_chr7$Fis)

ggsave("img/Fig_Fis.png",width = 6, height = 5)


