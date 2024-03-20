library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

here::i_am("Rproj_slccolonpaper/Final_Figures_Paper/Figure_7_DSS_TNBS.R")

# DSS Body Weight Curve
pct_weight <- read.csv(here("SLC_DSS/PCT_Body_Weight.csv"))
pct_weight_long <- pivot_longer(pct_weight, 
                                cols = starts_with("D"), 
                                names_to = "Day", 
                                values_to = "Score")

pct_weight_long$Day <- as.integer(stringr::str_extract(pct_weight_long$Day, "\\d+"))
pct_weight_long$Day <- factor(pct_weight_long$Day)
pct_weight_long$Score <- gsub("%","",pct_weight_long$Score)
pct_weight_long$Score <- as.numeric(pct_weight_long$Score)
pct_weight_long$Genotype <- factor(pct_weight_long$Genotype, levels=c("WT", "HET", "MUT"))
pct_weight_long <- remove_missing(pct_weight_long)

# Create plot
make_longitudinal_graph<- function(dataframe,ylab,xlab) {
  data_long<- as.data.frame(dataframe)
  df_summary <- data_long %>% 
    group_by(Genotype, Day) %>%
    summarise(mean = mean(Score), 
              sd = sd(Score),
              se = sd / sqrt(n()))
  
  ggplot(df_summary, aes(x = Day, y = mean, group = Genotype, color = Genotype)) +
    geom_line(size=2) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
    theme_cowplot(20) + 
    ggtitle(ylab) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
    theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")
}
percent_weight_dss <- make_longitudinal_graph(pct_weight_long,"DSS Body Weight (% Baseline)", "Day of DSS")

dss_summary <- pct_weight_long %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

dss_summary %>% kable

## TNBS Body Weight Curve --
pct_weight <- read.csv("TNBS/PCT_Body_Weight.csv")

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

percent_weight_tnbs<- make_longitudinal_graph(pct_weight_long,"TNBS Body Weight (% Baseline)", "Day of TNBS")

tnbs_summary <- pct_weight_long %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

tnbs_summary %>% kable


# Read / clean data
histology <- read.csv("SLC_DSS/Histology.csv")
histology <- remove_missing(histology)
histology$Genotype <- factor(histology$Genotype, levels=c("WT","HET","MUT"))

# Create plot
dss_histo_plot <- histology %>%
  ggplot( aes(x=Genotype, y=Score, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("DSS Histology Score") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

# Read / clean data
colon_length <- read.csv("TNBS/Colon_Length.csv")
colon_length <- remove_missing(colon_length)
colon_length$Genotype <- factor(colon_length$Genotype, levels=c("WT","HET","MUT"))

# Create plot
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

# Read / clean data
colon_length <- read.csv("SLC_DSS/Colon_and_Spleen.csv")
colon_length <- remove_missing(colon_length)
colon_length$Genotype <- factor(colon_length$Genotype, levels=c("WT","HET","MUT"))

# Create plot
dss_colon_plot <- colon_length %>%
  ggplot( aes(x=Genotype, y=Colon.Length, fill=Genotype)) +
  geom_boxplot(alpha=0.6) +
  scale_fill_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  geom_beeswarm(size=3)+
  theme_cowplot(20)+
  ggtitle("DSS Colon Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

# Plot arranged figures
plot_grid(dss_histo_plot, dss_colon_plot,
          percent_weight_tnbs, tnbs_colon_plot,
          nrow=2, label_size = 20,
          labels=c("A","B","C","D"))

