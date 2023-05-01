library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(tidyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Overlapping_Genera.R")

## Luminal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
data$Taxon <- data$feature
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub("\\..*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub("\\..*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub("\\..*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)

data$annotation <- gsub("\\.E","E",data$Genus)
data$annotation <- gsub("\\.","_",data$annotation)
data$annotation <- gsub("__","_",data$annotation)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

trios_l6_mut_significant<- data %>% filter(value=="MUT")
trios_l6_het_significant<- data %>% filter(value=="HET")

# Long term
data<-read.table("Long_Term/differential_taxa/L6_Luminal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
data$Taxon <- data$feature
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub("\\..*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub("\\..*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub("\\..*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)

data$annotation <- gsub("\\.E","E",data$Genus)
data$annotation <- gsub("\\.","_",data$annotation)
data$annotation <- gsub("__","_",data$annotation)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

long_term_l6_mut_significant <- data %>% filter(value=="MUT")
long_term_l6_het_significant <- data %>% filter(value=="HET")


# Baseline
data<-read.table("Baseline/differential_taxa/JAX_L6_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
data$Taxon <- data$feature
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub("\\..*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub("\\..*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub("\\..*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)

data$annotation <- gsub("\\.E","E",data$Genus)
data$annotation <- gsub("\\.","_",data$annotation)
data$annotation <- gsub("__","_",data$annotation)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

baseline_l6_mut_significant <- data %>% filter(value=="MUT")
baseline_l6_het_significant <- data %>% filter(value=="HET")

## LUMINAL HEATMAP ACROSS DATASETS --
baseline_l6_mut_significant$Dataset = c("Baseline_JAX")
long_term_l6_mut_significant$Dataset = c("Long-Term")
trios_l6_mut_significant$Dataset = c("Trios")

# number of intersecting features
v1 <- intersect(baseline_l6_mut_significant$feature, long_term_l6_mut_significant$feature) #32 features
v2 <- intersect(baseline_l6_mut_significant$feature, trios_l6_mut_significant$feature) #14 features
v3 <- intersect(trios_l6_mut_significant$feature, long_term_l6_mut_significant$feature) # 1 feature

v4 <- c(v1, v2,v3)
target <- unique(v4)

subset_baseline_l6_mut_significant <- baseline_l6_mut_significant[baseline_l6_mut_significant$feature %in% target, ]
subset_long_term_l6_mut_significant <-  long_term_l6_mut_significant[long_term_l6_mut_significant$feature %in% target, ]
subset_trios_l6_mut_significant <-  trios_l6_mut_significant[trios_l6_mut_significant$feature %in% target, ]

overlapping_l6_mut_luminal <- rbind(subset_baseline_l6_mut_significant, 
                                    subset_long_term_l6_mut_significant,
                                    subset_trios_l6_mut_significant)
df <- overlapping_l6_mut_luminal
# Create the heatmap
ggplot(df, aes(x = Dataset, y = annotation, fill = coef)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = "Dataset", y = "Feature") +
  theme_bw() +
  theme_cowplot(16)+
  ggtitle("Luminal Enriched MUT (yellow) vs Enriched WT (purple): 2 out of 3 datasets")+
  geom_text(aes(x = 0.5, y = annotation, label = annotation, color = Phylum), 
            size = 3, hjust = 1, fontface = "bold")# +

## Mucosal --
# Trios
data<-read.table("Trios/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
taxonomy <- read.delim("Long_Term/starting_files/final_taxonomy.tsv")
data$Taxon <- data$feature
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub("\\..*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub("\\..*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub("\\..*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)

data$annotation <- gsub("\\.E","E",data$Genus)
data$annotation <- gsub("\\.","_",data$annotation)
data$annotation <- gsub("__","_",data$annotation)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

trios_l6_mut_significant<- data %>% filter(value=="MUT")
trios_l6_het_significant<- data %>% filter(value=="HET")

# Long term
data<-read.table("Long_Term/differential_taxa/L6_Mucosal_Colon_L6_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")
taxonomy <- read.delim("Long_Term/starting_files/final_taxonomy.tsv")
data$Taxon <- data$feature
data$Phylum <- gsub(".*p__","",data$Taxon)
data$Phylum <- gsub("\\..*","",data$Phylum)
data$Family<- gsub(".*f__","",data$Taxon)
data$Family <-  gsub("\\..*","",data$Family)
data$Order<- gsub(".*o__","",data$Taxon)
data$Order <-  gsub("\\..*","",data$Order)
data$Genus<- gsub(".*g__","",data$Taxon)

data$annotation <- gsub("\\.E","E",data$Genus)
data$annotation <- gsub("\\.","_",data$annotation)
data$annotation <- gsub("__","_",data$annotation)
#data$Genus <- gsub("\\..*","",data$Genus)
data <- data %>% mutate(annotation = ifelse(data$Genus=="", paste0(data$Family," (f)"), data$annotation))
data <- data %>% mutate(annotation = ifelse(data$Family=="", paste(data$Order,"(o)"), data$annotation))

long_term_l6_mut_significant <- data %>% filter(value=="MUT")
long_term_l6_het_significant <- data %>% filter(value=="HET")

## MUCOSAL HEATMAP ACROSS DATASETS --
long_term_l6_mut_significant$Dataset = c("Long-Term")
trios_l6_mut_significant$Dataset = c("Trios")

# number of intersecting features
v3 <- intersect(trios_l6_mut_significant$feature, long_term_l6_mut_significant$feature) # 1 feature

target <-v3

subset_long_term_l6_mut_significant <-  long_term_l6_mut_significant[long_term_l6_mut_significant$feature %in% target, ]
subset_trios_l6_mut_significant <-  trios_l6_mut_significant[trios_l6_mut_significant$feature %in% target, ]

overlapping_l6_mut_mucosal <- rbind(subset_long_term_l6_mut_significant,
                                     subset_trios_l6_mut_significant)
df <- overlapping_l6_mut_mucosal
# Create the heatmap
ggplot(df, aes(x = Dataset, y = annotation, fill = coef)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = "Dataset", y = "Feature") +
  theme_bw() +
  theme_cowplot(16)+
  ggtitle("Mucosal Enriched MUT (yellow) vs Enriched WT (purple): shared Genera")+
  geom_text(aes(x = 0.5, y = annotation, label = annotation, color = Phylum), 
            size = 3, hjust = 1, fontface = "bold")# +
