library(ggplot2)
library(dplyr)
library(cowplot)
library(nlme)
library(tidyr)
library(knitr)

setwd("../SLC_DSS/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/SLC_DSS_phenotype.R")

data <- read.csv("Stool_Phenotype.csv",header=TRUE)

data_long <- pivot_longer(data, 
                          cols = starts_with("D"), 
                          names_to = "Day", 
                          values_to = "Score")

data_long$Day <- as.integer(stringr::str_extract(data_long$Day, "\\d+"))
data_long$Genotype <- factor(data_long$Genotype, levels=c("WT", "HET", "MUT"))
data_long <- remove_missing(data_long)

## All together -- 
# Calculate the mean and standard error for each group

make_longitudinal_graph <- function(filterby, ylab) {
df_summary <- data_long %>% filter(Phenotype==filterby) %>%
  group_by(Genotype, Day) %>%
  summarise(mean = mean(Score), 
            sd = sd(Score),
            se = sd / sqrt(n()))

# Plot the graph with error bars
ggplot(df_summary, aes(x = Day, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size=2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(x = "Day of DSS", y = ylab) +
  scale_color_manual(values = c("WT" ="black","MUT" = "red", "HET" = "blue")) +
  theme_cowplot(20) + 
  ggtitle(ylab) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(legend.position = "top",legend.title = element_text(hjust = 0.5), legend.justification = "center")

} 

clin_score <- make_longitudinal_graph("Clinical_Score", "Clinical Score")
stool_consist <- make_longitudinal_graph("Stool_Consistency", "Stool Consistency")
occult_score <- make_longitudinal_graph("Occult_Score", "Fecal Occult")

plot_grid(clin_score, stool_consist, occult_score, labels=c("C", "D", "E"), ncol=3, label_size = 22)


### Stats ---
data_long <- pivot_longer(data, 
                          cols = starts_with("D"), 
                          names_to = "Day", 
                          values_to = "Score")
clinical <- data_long %>% filter(Phenotype=="Clinical_Score")
occult <- data_long %>% filter(Phenotype=="Occult_Score")
consist <- data_long %>% filter(Phenotype=="Stool_Consistency")


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

## Splitting the dataframe by Day - Occult -
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
