
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(ggpubr)

## Environment --

here::i_am("Rproj_slccolonpaper/Trios_Alpha_Diversity.R")

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
otus<- readr::read_delim(here("Trios/alpha_diversity/alpha_min_10000_Luminal_Colon_notax_Trios_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv"))
row.names(otus) <- otus$...1
shannon<-readr::read_delim(here("Trios/alpha_diversity/alpha_min_10000_Luminal_Colon_notax_Trios_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv"))
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim(here("Trios/alpha_diversity/alpha_min_10000_Luminal_Colon_notax_Trios_ASV_table_Silva_v138_1/chao1_dir/alpha-diversity.tsv"))
row.names(chao1) <- chao1$...1

data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- readr::read_delim(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"))
data_meta <- merge(data,metadata, by="SampleID")
data_meta$Subset <- dplyr::recode(data_meta$Subset,Luminal_Colon = "Colon Lumen")

## Colon mucosa --
otus<- readr::read_delim(here("Trios/alpha_diversity/alpha_Mucosal_Colon_notax_Trios_ASV_table_Silva_v138_1/otus_dir/alpha-diversity.tsv"))
row.names(otus) <- otus$...1
shannon<-readr::read_delim(here("Trios/alpha_diversity/alpha_Mucosal_Colon_notax_Trios_ASV_table_Silva_v138_1/shannon_dir/alpha-diversity.tsv"))
row.names(shannon) <- shannon$...1
chao1<-readr::read_delim(here("Trios/alpha_diversity/alpha_Mucosal_Colon_notax_Trios_ASV_table_Silva_v138_1/chao1_dir/alpha-diversity.tsv"))
row.names(chao1) <- chao1$...1

data<- merge(otus,shannon, by="...1")
data<- merge(data,chao1, by="...1")
data$SampleID <- data$...1

metadata<- readr::read_delim(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"))
mucosal_data_meta <- merge(data,metadata, by="SampleID")
mucosal_data_meta$Subset <- dplyr::recode(data_meta$Subset,Mucosal_Colon = "Colon Mucosa")

### Make and store plots ---
adiv_trios_shannon<- generate_adiv_plots(data_meta, Genotype, shannon_entropy, 2, 7) +
  facet_wrap(~Subset)+
  labs(x="",y="shannon")+
  ggtitle("2 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_trios_shannon

adiv_trios_otus<- generate_adiv_plots(data_meta, Genotype, observed_features, 0, 300) +
  facet_wrap(~Subset)+
  labs(x="",y="# ASVs")+
  ggtitle("2 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
adiv_trios_otus

mc_adiv_trios_shannon<- generate_adiv_plots(mucosal_data_meta, Genotype, shannon_entropy, 2, 7) +
  facet_wrap(~Subset)+
  labs(x="",y="shannon")+
  ggtitle("2 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
mc_adiv_trios_shannon

mc_adiv_trios_otus<- generate_adiv_plots(mucosal_data_meta, Genotype, observed_features, 0, 300) +
  facet_wrap(~Subset)+
  labs(x="",y="# ASVs")+
  ggtitle("2 month-old")+
  theme(plot.title = element_text(hjust = 0.5))
mc_adiv_trios_otus


plot_grid(adiv_trios_shannon, adiv_trios_otus, adiv_trios_chao1, labels=c("A","B","C"), nrow=1)

### Alpha Diversity Stats ---
data_meta$Genotype <-factor(data_meta$Genotype, levels=c("WT", "HET","MUT"))
output <- lme(fixed= shannon_entropy ~ Sequencing_Run + Sex + Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sequencing_Run +Sex+ Site+ Genotype, random = ~1|MouseID, data=data_meta)
summary(output)

mucosal_data_meta$Genotype <-factor(mucosal_data_meta$Genotype, levels=c("WT", "HET","MUT"))
output <- lme(fixed= shannon_entropy ~ Sequencing_Run + Sex + Site+ Genotype, random = ~1|MouseID, data=mucosal_data_meta)
summary(output)
output <- lme(fixed= observed_features ~ Sequencing_Run +Sex+ Site+ Genotype, random = ~1|MouseID, data=mucosal_data_meta)
summary(output)

