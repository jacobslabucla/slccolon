library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/SLC_Spontaneous_Phenotype.R")

mucin <- read.csv("Spontaneous/Mucin.csv")
histology <- read.csv("Spontaneous/Histology.csv")

histology$Genotype <- factor(histology$Genotype, levels=c("WT","HET","MUT"))

histo_plot <- histology %>%
  ggplot( aes(x=Genotype, y=Score, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("Histology Score") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

mucin$Genotype <- factor(mucin$Genotype, levels=c("WT","HET","MUT"))

mucin_plot <- mucin %>%
  ggplot( aes(x=Genotype, y=Mucin_Thickness, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("Mucin Thickness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

plot_grid(histo_plot, NULL, 
          mucin_plot, NULL, 
          labels=c("A", "B", "C", "D"),
          label_size=20)
