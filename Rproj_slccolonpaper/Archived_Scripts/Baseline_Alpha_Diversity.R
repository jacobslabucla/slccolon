library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(ggpubr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolon/")

## Environment --
here::i_am("Rproj_slccolonpaper/Baseline_Alpha_Diversity.R")

otus <- readr::read_delim(here("Baseline/alpha_diversity/alpha_min_10000_Baseline_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv"))
row.names(otus) <- otus$...1
shannon<-readr::read_delim(here("Baseline/alpha_diversity/alpha_min_10000_Baseline_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv"))
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim(here("Baseline/alpha_diversity/alpha_min_10000_Baseline_ASV_table_Silva_v138_1/chao1_dir/alpha-diversity.tsv"))
row.names(chao1) <- chao1$...1

data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- readr::read_csv(here("Baseline/starting_files/Baseline_Metadata.csv"))
metadata$SampleID <- metadata$...1
data_meta <- merge(data,metadata, by="SampleID")
data_meta <- data_meta %>% filter(Line!="Q22")

### Function for plotting alpha diversity ---
generate_adiv_plots <- function(input_data, X, Y, min, max){
  #read in files
  data <- as.data.frame(input_data)
  
  #declare order of variables
  data$Genotype <- factor(data$Genotype, levels=c("WT", "HET","MUT"))
  #graph plot
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill=Genotype)) + 
    geom_boxplot(alpha=0.25)+
    scale_fill_viridis_d()+
    #geom_line(aes(group = MouseID,color=Genotype),size=1)+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(16) +
    #ylim(min,max) +
    theme(legend.position = "none")
  
}

### Make and store plots ---
compare <-c(c("WT","HET"), c("WT","MUT"))

baseline_shannon<- generate_adiv_plots(data_meta, Genotype, shannon, 0, 7) +
  stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  ggtitle("3-4 months")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
baseline_shannon

baseline_otus<- generate_adiv_plots(data_meta, Genotype, observed_otus, 0, 300) +
  ylab("# ASVs") +
  xlab("")+
  ggtitle("3-4 months")+
  theme(plot.title = element_text(hjust = 0.5))
baseline_otus

