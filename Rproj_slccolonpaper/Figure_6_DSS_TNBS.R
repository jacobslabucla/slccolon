library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Figure_6_DSS_TNBS.R")

plot_grid(percent_weight, NULL, dss_histo_plot, dss_colon_plot,
          NULL, NULL, percent_weight_tnbs, tnbs_colon_plot,
          nrow=2, label_size = 20,
          labels=c("A", "B","C","D","","","E","F"))


## TNBS colon length -- 
colon_length <- read.csv("TNBS/Colon_Length.csv")
colon_length <- remove_missing(colon_length)
colon_length$Genotype <- factor(colon_length$Genotype, levels=c("WT","HET","MUT"))

tnbs_colon_plot <- colon_length %>%
  ggplot( aes(x=Genotype, y=Colon_Length, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("TNBS Colon Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")
