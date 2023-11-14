library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(tidyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Overlapping_ECs.R")

## LUMINAL EC COMPARISONS --

# Trios
data<-read.table("Trios/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Trios/differential_EC/annotated_EC.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
annotation$feature
annotation$enzyme_class <- sapply(annotation$feature, assign_letter)
data <- merge(data,annotation, by="feature")

trios_ec_mut_significant_LC <- data %>% filter(value=="MUT")
trios_ec_het_significant_LC <- data %>% filter(value=="HET")

# Long term
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Long_Term/differential_EC/annotated_ec.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
annotation$feature
annotation$enzyme_class <- sapply(annotation$feature, assign_letter)
data <- merge(data,annotation, by="feature")

long_term_ec_mut_significant <- data %>% filter(value=="MUT")
long_term_ec_het_significant <- data %>% filter(value=="HET")


# Baseline
data<-read.table("Baseline/Baseline_JAX_EC_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Baseline/differential_EC/annotated_ec.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
annotation$feature
annotation$enzyme_class <- sapply(annotation$feature, assign_letter)
data <- merge(data,annotation, by="feature")

baseline_ec_mut_significant <- data %>% filter(value=="MUT")
baseline_ec_het_significant <- data %>% filter(value=="HET")

## LUMINAL HEATMAP ACROSS DATASETS --
baseline_ec_mut_significant$Dataset = c("Baseline_JAX")
long_term_ec_mut_significant$Dataset = c("Long-Term")
trios_ec_mut_significant_LC$Dataset = c("Trios")

# number of intersecting features
v1 <- intersect(baseline_ec_mut_significant$feature, long_term_ec_mut_significant$feature) #32 features
v2 <- intersect(baseline_ec_mut_significant$feature, trios_ec_mut_significant_LC$feature) #14 features
v3 <- intersect(trios_ec_mut_significant_LC$feature, long_term_ec_mut_significant$feature) # 1 feature

v4 <- c(v1, v2,v3)
target <- unique(v4)

subset_baseline_ec_mut_significant <- baseline_ec_mut_significant[baseline_ec_mut_significant$feature %in% target, ]
subset_long_term_ec_mut_significant <-  long_term_ec_mut_significant[long_term_ec_mut_significant$feature %in% target, ]
subset_trios_ec_mut_significant <-  trios_ec_mut_significant_LC[trios_ec_mut_significant_LC$feature %in% target, ]

overlapping_ec_mut_luminal <- rbind(subset_baseline_ec_mut_significant, 
                                    subset_long_term_ec_mut_significant,
                                    subset_trios_ec_mut_significant)
df <- overlapping_ec_mut_luminal
# Create the heatmap
ggplot(df, aes(x = Dataset, y = description, fill = coef)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = "Dataset", y = "Feature") +
  theme_bw() +
  theme_cowplot(16)+
  ggtitle("Luminal Enriched MUT (yellow) vs Enriched WT (purple): 2 out of 3 datasets")+
  geom_text(aes(x = 0.5, y = description, label = description, color = enzyme_class), 
            size = 3, hjust = 1, fontface = "bold")# +

## MUCOSAL EC COMPARISONS --

# Trios
data<-read.table("Trios/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Trios/differential_EC/annotated_EC.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
annotation$feature
annotation$enzyme_class <- sapply(annotation$feature, assign_letter)
data <- merge(data,annotation, by="feature")

trios_ec_mut_significant_LC <- data %>% filter(value=="MUT")
trios_ec_het_significant_LC <- data %>% filter(value=="HET")

# Long term
data<-read.table("Long_Term/differential_EC/PICRUST2_EC_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
annotation <- read.delim("Long_Term/differential_EC/annotated_ec.tsv", row.names=1)
annotation$feature <- row.names(annotation)
annotation <- annotation %>% select(c("feature","description"))
annotation$feature <- gsub(":", ".", annotation$feature)
annotation$feature
annotation$enzyme_class <- sapply(annotation$feature, assign_letter)
data <- merge(data,annotation, by="feature")

long_term_ec_mut_significant <- data %>% filter(value=="MUT")
long_term_ec_het_significant <- data %>% filter(value=="HET")

## MUCOSAL HEATMAP ACROSS DATASETS --
long_term_ec_mut_significant$Dataset = c("Long-Term")
trios_ec_mut_significant_LC$Dataset = c("Trios")

# number of intersecting features
v3 <- intersect(trios_ec_mut_significant_LC$feature, long_term_ec_mut_significant$feature) # 7 feature

target <- unique(v3)

subset_long_term_ec_mut_significant <-  long_term_ec_mut_significant[long_term_ec_mut_significant$feature %in% target, ]
subset_trios_ec_mut_significant <-  trios_ec_mut_significant_LC[trios_ec_mut_significant_LC$feature %in% target, ]

overlapping_ec_mut_mucosal <- rbind(subset_long_term_ec_mut_significant,
                                    subset_trios_ec_mut_significant)
df <- overlapping_ec_mut_mucosal

# Create the heatmap
ggplot(df, aes(x = Dataset, y = description, fill = coef)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = "Dataset", y = "Feature") +
  theme_bw() +
  theme_cowplot(16)+
  ggtitle("Mucosal Enriched MUT (yellow) vs Enriched WT (purple): shared in 2 datasets")+
  geom_text(aes(x = 0.5, y = description, label = description, color = enzyme_class), 
            size = 3, hjust = 1, fontface = "bold")# +

