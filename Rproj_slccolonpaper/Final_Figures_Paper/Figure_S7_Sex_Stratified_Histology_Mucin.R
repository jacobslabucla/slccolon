library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)
library(here)

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_S7_Sex_Stratified_Histology_Mucin.R")

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
  theme(legend.position = "none",legend.title = element_text(hjust = 0.5), legend.justification = "center")+
  facet_wrap(~Sex)

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
  theme(legend.position = "none",legend.title = element_text(hjust = 0.5), legend.justification = "center")+
  facet_wrap(~Sex)

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
  theme(legend.position = "none",legend.title = element_text(hjust = 0.5), legend.justification = "center")+
  facet_wrap(~Sex)

plot_grid(histo_plot,histo_old_plot, mucin_plot, NULL,
          labels=c("A","B","C",""),
          nrow=2,
          label_size=12)


## Statistics -- 

# 5 month
wilcox.test(Score~Genotype, histology)

f_df_mut <- histology %>% filter(Sex=="Female")
wilcox.test(Score~Genotype, f_df_mut)

m_df_mut <- histology %>% filter(Sex=="Male")
wilcox.test(Score~Genotype, m_df_mut)

# 10 month

df_mut <- histo_old %>% filter(Genotype!="HET")
wilcox.test(Score~Genotype, df_mut)

f_df_mut <- histo_old %>% filter(Genotype!="HET" & Sex=="Female")
wilcox.test(Score~Genotype, f_df_mut)

m_df_mut <- histo_old %>% filter(Genotype!="HET" & Sex=="Male")
wilcox.test(Score~Genotype, m_df_mut)

# Mucin 
df_het <- mucin %>% filter(Genotype!="MUT")
wilcox.test(Mucin_Thickness~Genotype, df_het)

f_df_het <- mucin %>% filter(Genotype!="MUT" & Sex =="Female")
wilcox.test(Mucin_Thickness~Genotype, f_df_het)

m_df_het <- mucin %>% filter(Genotype!="MUT" & Sex =="Male")
wilcox.test(Mucin_Thickness~Genotype, m_df_het)

df_mut <- mucin %>% filter(Genotype!="HET")
wilcox.test(Mucin_Thickness~Genotype,df_mut)

f_df_mut <- mucin %>% filter(Genotype!="HET" & Sex=="Female")
wilcox.test(Mucin_Thickness~Genotype,f_df_mut)

m_df_mut <- mucin %>% filter(Genotype!="HET" & Sex=="Male")
wilcox.test(Mucin_Thickness~Genotype,m_df_mut)

df_mut <- mucin %>% filter(Genotype!="WT")
wilcox.test(Mucin_Thickness~Genotype,df_mut)

f_df_mut <- mucin %>% filter(Genotype!="WT")
wilcox.test(Mucin_Thickness~Genotype,f_df_mut)

m_df_mut <- mucin %>% filter(Genotype!="WT")
wilcox.test(Mucin_Thickness~Genotype,m_df_mut)
