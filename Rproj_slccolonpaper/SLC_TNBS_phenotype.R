library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/SLC_TNBS_phenotype.R")

weight <- read.csv("TNBS/Weight.csv")
pct_weight <- read.csv("TNBS/PCT_Body_Weight.csv")
colon_length <- read.csv("TNBS/Colon_Length.csv")

## Graph Colon Length -- 
colon_length$Genotype <- factor(colon_length$Genotype,levels=c("WT","HET","MUT"))
colon_plot <- colon_length %>%
  ggplot( aes(x=Genotype, y=Colon_Length, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ylab("Colon Length (cm)")+
  ggtitle("TNBS Colon Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

## Wrangle longitudinal data into long formats --
weight <- weight %>% filter(Batch!="Three")
weight_long <- pivot_longer(weight, 
                            cols = starts_with("D"), 
                            names_to = "Day", 
                            values_to = "Score")

weight_long$Day <- as.integer(stringr::str_extract(weight_long$Day, "\\d+"))
weight_long$Day <- factor(weight_long$Day)
weight_long$Genotype <- factor(weight_long$Genotype, levels=c("WT", "HET", "MUT"))
weight_long <- remove_missing(weight_long)

pct_weight_long <- pct_weight %>% filter(Batch!="Three")
pct_weight_long <- pivot_longer(pct_weight, 
                                cols = starts_with("D"), 
                                names_to = "Day", 
                                values_to = "Score")

pct_weight_long$Day <- as.integer(stringr::str_extract(pct_weight_long$Day, "\\d+"))
pct_weight_long$Day <- factor(pct_weight_long$Day)
pct_weight_long$Score <- as.numeric(pct_weight_long$Score)*100
pct_weight_long$Score <- as.numeric(pct_weight_long$Score)
pct_weight_long$Genotype <- factor(pct_weight_long$Genotype, levels=c("WT", "HET", "MUT"))
pct_weight_long <- remove_missing(pct_weight_long)

## All together -- 
# Calculate the mean and standard error for each group
make_longitudinal_graph_TNBS<- function(dataframe,ylab) {
  data_long<- as.data.frame(dataframe)
  df_summary <- data_long %>% 
    group_by(Genotype, Day) %>%
    summarise(mean = mean(Score), 
              sd = sd(Score),
              se = sd / sqrt(n()))
  
  # Plot the graph with error bars
  ggplot(df_summary, aes(x = Day, y = mean, group = Genotype, color = Genotype)) +
    geom_line(size=2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(x = "Day of TNBS", y = ylab) +
    scale_color_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
    theme_cowplot(20) + 
    ggtitle(ylab) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
    theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")
  
} 
make_longitudinal_graph_TNBS_batch <- function(dataframe,ylab) {
  data_long<- as.data.frame(dataframe)
  df_summary <- data_long %>% 
    group_by(Batch, Genotype, Day) %>%
    summarise(mean = mean(Score), 
              sd = sd(Score),
              se = sd / sqrt(n()))
  
  # Plot the graph with error bars
  ggplot(df_summary, aes(x = Day, y = mean, group = Genotype, color = Genotype)) +
    geom_line(size=2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(x = "Day of TNBS", y = ylab) +
    scale_color_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
    theme_cowplot(20) + 
    ggtitle(ylab) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
    theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")
  
} 

weight_raw_tnbs <- make_longitudinal_graph_TNBS(weight_long, "Body Weight (g)")
percent_weight_tnbs<- make_longitudinal_graph_TNBS(pct_weight_long,"TNBS Body Weight (% Baseline)")

plot_grid(weight_raw_tnbs, percent_weight_tnbs,labels=c("A","B"), label_size = 22)

weight_raw_tnbs_batch <- make_longitudinal_graph_TNBS_batch(weight_long, "Body Weight (g)") + facet_grid(~Batch)
percent_weight_tnbs_batch<- make_longitudinal_graph_TNBS_batch(pct_weight_long,"Body Weight (% Baseline)") + facet_grid(~Batch)

### Stats ---

# Pairwise compairison -

wilcox_test_to_markdown(pct_weight_long, "TNBS/PCT_Weight_Wilcoxon.md")

output <- lme(fixed= Score ~ Sex + Day*Genotype, random = ~1|MouseID, data=pct_weight_long)
write_lme_summary_to_md(output, "TNBS/DayANDGenotype_PCT_Weight_LMEM_summary.md")

output <- lme(fixed= Score ~ Sex + Day + Genotype, random = ~1|MouseID, data=pct_weight_long)
write_lme_summary_to_md(output, "TNBS/Genotype_PCT_Weight_LMEM_summary.md")

## Colon Length Stats --
df_het <- colon_length %>% filter(Genotype!="MUT")
wilcox.test(Colon_Length~Genotype, df_het)

df_mut <- colon_length %>% filter(Genotype!="HET")
wilcox.test(Colon_Length~Genotype,df_mut)

