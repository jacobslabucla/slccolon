library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)
library(here)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
#setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") 
here::i_am("Rproj_slccolonpaper/SLC_Spontaneous_phenotype.R")

mucin <- read.csv("Spontaneous/Mucin.csv")
histology <- read.csv("Spontaneous/Histology_5month.csv")
histo_old <- read.csv("Spontaneous/Histology_10month.csv")

histology$Genotype <- factor(histology$Genotype, levels=c("WT","HET","MUT"))

histo_plot <- histology %>%
  ggplot( aes(x=Genotype, y=Score, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=1)+
  theme_cowplot(12)+
  ggtitle("Histology Score") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

histo_old$Genotype <- factor(histo_old$Genotype, levels=c("WT","HET","MUT"))
histo_old <- histo_old %>% filter(Tg=="Negative")
histo_old_plot <- histo_old %>%
  ggplot( aes(x=Genotype, y=Score, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=1)+
  theme_cowplot(12)+
  ggtitle("Histology Score") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

mucin$Genotype <- factor(mucin$Genotype, levels=c("WT","HET","MUT"))

mucin_plot <- mucin %>%
  ggplot( aes(x=Genotype, y=Mucin_Thickness, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=1)+
  theme_cowplot(12)+
  ggtitle("Mucin Thickness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

mucin_plot <- mucin %>%
  ggplot( aes(x=Genotype, y=Mucin_Thickness, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("Mucin Thickness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")+
  facet_wrap(~Sex)

plot_grid(NULL, NULL, histo_plot, 
          NULL, NULL, NULL, 
          labels=c("A", "", "B","C", "","D"),
          nrow=2,
          label_size=12)

plot_grid(mucin_plot, NULL, 
          nrow=2, 
          labels=c("A","B"),
          label_size=12)


## Statistics -- 
df_het <- mucin %>% filter(Genotype!="MUT")
wilcox.test(Mucin_Thickness~Genotype, df_het)

df_mut <- mucin %>% filter(Genotype!="HET")
wilcox.test(Mucin_Thickness~Genotype,df_mut)

df_mut <- mucin %>% filter(Genotype!="WT")
wilcox.test(Mucin_Thickness~Genotype,df_mut)

df_mut <- histo_old %>% filter(Genotype!="HET")
wilcox.test(Score~Genotype, df_mut)

