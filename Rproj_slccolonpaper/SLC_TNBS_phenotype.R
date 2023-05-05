library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)
library(ggbeeswarm)

setwd("../") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/SLC_TNBS_phenotype.R")

weight <- read.csv("TNBS/Weight.csv")
pct_weight <- read.csv("TNBS/PCT_Body_Weight.csv")

## Wrangle longitudinal data into long formats --

weight_long <- pivot_longer(weight, 
                            cols = starts_with("D"), 
                            names_to = "Day", 
                            values_to = "Score")

weight_long$Day <- as.integer(stringr::str_extract(weight_long$Day, "\\d+"))
weight_long$Day <- factor(weight_long$Day)
weight_long$Genotype <- factor(weight_long$Genotype, levels=c("WT", "HET", "MUT"))
weight_long <- remove_missing(weight_long)

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

make_longitudinal_graph_TNBS <- function(dataframe,ylab) {
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

weight_raw_tnbs <- make_longitudinal_graph_TNBS(weight_long, "Body Weight (g)")
percent_weight_tnbs<- make_longitudinal_graph_TNBS(pct_weight_long,"Body Weight (% Baseline)")

plot_grid(weight_raw_tnbs, percent_weight_tnbs,labels=c("A","B"), label_size = 22)


### Stats ---

weight_long <- pivot_longer(weight, 
                            cols = starts_with("D"), 
                            names_to = "Day", 
                            values_to = "Score")

pct_weight_long <- pivot_longer(pct_weight, 
                                cols = starts_with("D"), 
                                names_to = "Day", 
                                values_to = "Score")
pct_weight_long$Score <- gsub("%","",pct_weight_long$Score)
pct_weight_long$Score <- as.numeric(pct_weight_long$Score)
## Pairwise compairison --

## Splitting the dataframe by Day - Clinical -
df_list <- split(clinical, clinical$Day)

# Defining the function to perform Wilcoxon rank sum test
wilcox_test <- function(df) {
  wt_mut <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "MUT"),])
  wt_het <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "HET"),])
  return(list(wt_mut, wt_het))
}

# Applying the function to each split dataframe
results <- lapply(df_list, wilcox_test)
results[1]
# Combining the results into a data frame
results_df <- do.call(rbind, lapply(seq_along(df_list), function(i) {
  day <- names(df_list)[i]
  res <- results[[i]]
  data.frame(Day = rep(day, 4),
             Genotype = rep(c("WT-MUT", "MUT-WT", "WT-HET", "HET-WT"), each = 1),
             W = c(res[[1]]$statistic, res[[1]]$statistic, res[[2]]$statistic, res[[2]]$statistic),
             p.value = c(res[[1]]$p.value, res[[1]]$p.value, res[[2]]$p.value, res[[2]]$p.value))
}))

# Convert the data frame to a markdown table
md_table <- knitr::kable(results_df, format = "markdown", align = c("c", "c", "c", "c"))

# Write the markdown table to a file
writeLines(md_table, "DSS_Clinical_Score_wilcoxon_results.md")

## Splitting the dataframe by Day - Occult --
df_list <- split(occult, occult$Day)

# Defining the function to perform Wilcoxon rank sum test
wilcox_test <- function(df) {
  wt_mut <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "MUT"),])
  wt_het <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "HET"),])
  return(list(wt_mut, wt_het))
}

# Applying the function to each split dataframe
results <- lapply(df_list, wilcox_test)

# Combining the results into a data frame
results_df <- do.call(rbind, lapply(seq_along(df_list), function(i) {
  day <- names(df_list)[i]
  res <- results[[i]]
  data.frame(Day = rep(day, 4),
             Genotype = rep(c("WT-MUT", "MUT-WT", "WT-HET", "HET-WT"), each = 1),
             W = c(res[[1]]$statistic, res[[1]]$statistic, res[[2]]$statistic, res[[2]]$statistic),
             p.value = c(res[[1]]$p.value, res[[1]]$p.value, res[[2]]$p.value, res[[2]]$p.value))
}))

# Convert the data frame to a markdown table
md_table <- knitr::kable(results_df, format = "markdown", align = c("c", "c", "c", "c"))

# Write the markdown table to a file
writeLines(md_table, "DSS_Occult_Score_wilcoxon_results.md")

## Splitting the dataframe by Day - Stool -
df_list <- split(consist, consist$Day)

# Defining the function to perform Wilcoxon rank sum test
wilcox_test <- function(df) {
  wt_mut <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "MUT"),])
  wt_het <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "HET"),])
  return(list(wt_mut, wt_het))
}

# Applying the function to each split dataframe
results <- lapply(df_list, wilcox_test)

# Combining the results into a data frame
results_df <- do.call(rbind, lapply(seq_along(df_list), function(i) {
  day <- names(df_list)[i]
  res <- results[[i]]
  data.frame(Day = rep(day, 4),
             Genotype = rep(c("WT-MUT", "MUT-WT", "WT-HET", "HET-WT"), each = 1),
             W = c(res[[1]]$statistic, res[[1]]$statistic, res[[2]]$statistic, res[[2]]$statistic),
             p.value = c(res[[1]]$p.value, res[[1]]$p.value, res[[2]]$p.value, res[[2]]$p.value))
}))

# Convert the data frame to a markdown table
md_table <- knitr::kable(results_df, format = "markdown", align = c("c", "c", "c", "c"))

# Write the markdown table to a file
writeLines(md_table, "DSS_Stool_Consistency_wilcoxon_results.md")

## Splitting the dataframe by Day - Clinical -
df_list <- split(weight_long, weight_long$Day)

# Defining the function to perform Wilcoxon rank sum test
wilcox_test <- function(df) {
  wt_mut <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "MUT"),])
  wt_het <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "HET"),])
  return(list(wt_mut, wt_het))
}

# Applying the function to each split dataframe
results <- lapply(df_list, wilcox_test)
results[1]
# Combining the results into a data frame
results_df <- do.call(rbind, lapply(seq_along(df_list), function(i) {
  day <- names(df_list)[i]
  res <- results[[i]]
  data.frame(Day = rep(day, 4),
             Genotype = rep(c("WT-MUT", "MUT-WT", "WT-HET", "HET-WT"), each = 1),
             W = c(res[[1]]$statistic, res[[1]]$statistic, res[[2]]$statistic, res[[2]]$statistic),
             p.value = c(res[[1]]$p.value, res[[1]]$p.value, res[[2]]$p.value, res[[2]]$p.value))
}))

# Convert the data frame to a markdown table
md_table <- knitr::kable(results_df, format = "markdown", align = c("c", "c", "c", "c"))

# Write the markdown table to a file
writeLines(md_table, "SLC_DSS/BW_raw_wilcoxon_results.md")

## Splitting the dataframe by Day - Clinical -
df_list <- split(pct_weight_long, pct_weight_long$Day)

# Defining the function to perform Wilcoxon rank sum test
wilcox_test <- function(df) {
  wt_mut <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "MUT"),])
  wt_het <- wilcox.test(Score ~ Genotype, data = df[df$Genotype %in% c("WT", "HET"),])
  return(list(wt_mut, wt_het))
}

# Applying the function to each split dataframe
results <- lapply(df_list, wilcox_test)
results[1]
# Combining the results into a data frame
results_df <- do.call(rbind, lapply(seq_along(df_list), function(i) {
  day <- names(df_list)[i]
  res <- results[[i]]
  data.frame(Day = rep(day, 4),
             Genotype = rep(c("WT-MUT", "MUT-WT", "WT-HET", "HET-WT"), each = 1),
             W = c(res[[1]]$statistic, res[[1]]$statistic, res[[2]]$statistic, res[[2]]$statistic),
             p.value = c(res[[1]]$p.value, res[[1]]$p.value, res[[2]]$p.value, res[[2]]$p.value))
}))

# Convert the data frame to a markdown table
md_table <- knitr::kable(results_df, format = "markdown", align = c("c", "c", "c", "c"))

# Write the markdown table to a file
writeLines(md_table, "SLC_DSS/PCT_BW_wilcoxon_results.md")

## Histo Stats --
df_het <- histology %>% filter(Genotype!="MUT")
wilcox.test(Score~Genotype, df_het)

df_mut <- histology %>% filter(Genotype!="HET")
wilcox.test(Score~Genotype,df_mut)

## LMEM -- 
output <- lme(fixed= FP_output ~ Sex + SLC_Genotype + ASO_Tg, random = ~1|MouseID, data=data_long)
summary(output)

pos <- data_long %>% filter(ASO_Tg=="Positive")
neg <- data_long %>% filter(ASO_Tg=="Negative")

output <- lme(fixed= FP_output ~ Sex +SLC_Genotype, random = ~1|MouseID, data=pos)
summary(output)
output <- lme(fixed= FP_output ~ Sex+SLC_Genotype, random = ~1|MouseID, data=neg)
summary(output)

pos_males <- data_long %>% filter(ASO_Tg=="Positive"& Sex=="Male")
neg_males <- data_long %>% filter(ASO_Tg=="Negative"& Sex=="Male")

output <- lme(fixed= FP_output ~ SLC_Genotype, random = ~1|MouseID, data=pos_males)
summary(output)
output <- lme(fixed= FP_output ~ SLC_Genotype, random = ~1|MouseID, data=neg_males)
summary(output)

pos_females <- data_long %>% filter(ASO_Tg=="Positive"& Sex=="Female")
neg_females <- data_long %>% filter(ASO_Tg=="Negative"& Sex=="Female")

output <- lme(fixed= FP_output ~ SLC_Genotype, random = ~1|MouseID, data=pos_females)
summary(output)
output <- lme(fixed= FP_output ~ SLC_Genotype, random = ~1|MouseID, data=neg_females)
summary(output)
