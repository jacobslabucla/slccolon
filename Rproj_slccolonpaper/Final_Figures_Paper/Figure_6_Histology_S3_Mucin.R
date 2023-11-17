library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)
library(here)

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_6_Histology_S3_Mucin.R")

mucin <- readr::read_csv(here("Spontaneous/Mucin.csv"))
histology <- readr::read_csv(here("Spontaneous/Histology_5month.csv"))
histo_old <- readr::read_csv(here("Spontaneous/Histology_10month.csv"))

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

plot_grid(histo_plot,histo_old_plot, 
          labels=c("B","D"),
          nrow=2,
          label_size=12)


## Statistics -- 
df_mut <- histo_old %>% filter(Genotype!="HET")
wilcox.test(Score~Genotype, df_mut)

## Statistics for mucin-- 
df_het <- mucin %>% filter(Genotype!="MUT")
wilcox.test(Mucin_Thickness~Genotype, df_het)
output <- lm(Mucin_Thickness~Sex+Genotype, df_het)
summary(output)

df_mut <- mucin %>% filter(Genotype!="HET")
wilcox.test(Mucin_Thickness~Genotype,df_mut)
output <- lm(Mucin_Thickness~Sex+Genotype, df_mut)
summary(output)

df_mut <- mucin %>% filter(Genotype!="WT")
wilcox.test(Mucin_Thickness~Genotype,df_mut)