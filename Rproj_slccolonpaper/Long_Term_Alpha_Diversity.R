
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(ggpubr)

## Environment --

here::i_am("Rproj_slccolonpaper/Baseline_Alpha_Diversity.R")

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

## Colon lumen --
otus<- readr::read_delim(here("Long_Term/alpha_diversity/alpha_Luminal_Colon_SLT_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv"))
row.names(otus) <- otus$...1
shannon<-readr::read_delim(here("Long_Term/alpha_diversity/alpha_Luminal_Colon_SLT_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv"))
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim(here("Long_Term/alpha_diversity/alpha_Luminal_Colon_SLT_ASV_table_Silva_v138_1/chao1_dir/alpha-diversity.tsv"))
row.names(chao1) <- chao1$...1


data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- readr::read_delim(here("Long_Term/starting_files/SLC_LT_metadata.tsv"))
data_meta <- merge(data,metadata, by="SampleID")
data_meta$Subset <- dplyr::recode(data_meta$Subset,Luminal_Colon = "Colon Lumen")

## Colon mucosa
otus<- readr::read_delim(here("Long_Term/alpha_diversity/alpha_Mucosal_Colon_SLT_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv"))
row.names(otus) <- otus$...1
shannon<-readr::read_delim(here("Long_Term/alpha_diversity/alpha_Mucosal_Colon_SLT_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv"))
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim(here("Long_Term/alpha_diversity/alpha_Mucosal_Colon_SLT_ASV_table_Silva_v138_1/chao1_dir/alpha-diversity.tsv"))
row.names(chao1) <- chao1$...1

data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- readr::read_delim(here("Long_Term/starting_files/SLC_LT_metadata.tsv"))
mucosal_data_meta <- merge(data,metadata, by="SampleID")
mucosal_data_meta$Subset <- dplyr::recode(data_meta$Subset,Mucosal_Colon = "Colon Mucosa")

### Make and store plots ---
compare <-c(c("WT","HET"), c("WT","MUT"))

adiv_slt_shannon<- generate_adiv_plots(data_meta, Genotype, shannon_entropy, 2, 7) +
  #stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  facet_wrap(~Subset)+
  labs(x="",y="shannon")+
  ggtitle("12 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_slt_shannon

adiv_slt_otus<- generate_adiv_plots(data_meta, Genotype, observed_features, 0, 300) +
  #stat_compare_means(comparisons = compare,method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)+
  facet_wrap(~Subset)+
  labs(x="",y="# ASVs")+
  ggtitle("12 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_slt_otus


mc_adiv_slt_shannon<- generate_adiv_plots(mucosal_data_meta, Genotype, shannon_entropy, 2, 7) +
  facet_wrap(~Subset)+
  labs(x="",y="shannon")+
  ggtitle("12 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
mc_adiv_slt_shannon

mc_adiv_slt_otus<- generate_adiv_plots(mucosal_data_meta, Genotype, observed_features, 0, 300) +
  facet_wrap(~Subset)+
  labs(x="",y="# ASVs")+
  ggtitle("12 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
mc_adiv_slt_otus

plot_grid(adiv_slt_shannon, adiv_slt_otus, mc_adiv_slt_shannon,mc_adiv_slt_otus, labels=c("A","B","C","D"), nrow=2)

### Alpha Diversity Stats ---
data_meta$Genotype <-factor(data_meta$Genotype, levels=c("WT", "HET","MUT"))
output <- lme(fixed= shannon_entropy ~ Sex + Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sex+ Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= chao1 ~ Sex+ Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)

mucosal_data_meta$Genotype <-factor(mucosal_data_meta$Genotype, levels=c("WT", "HET","MUT"))
output <- lme(fixed= shannon_entropy ~ Sex + Site+ Genotype, random = ~1|MouseID, data=mucosal_data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sex+ Site+ Genotype, random = ~1|MouseID, data=mucosal_data_meta)
summary(output)


