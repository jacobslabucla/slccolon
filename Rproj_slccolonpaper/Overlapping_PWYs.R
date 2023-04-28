library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/")

# function to load and process data
load_and_process_data <- function(filepath_to_significant_results, filepath_to_annotation) {
  # load the significant results
  data <- read.table(filepath_to_significant_results, header=TRUE)
  data <- data %>% filter(qval < 0.25)
  data <- data %>% filter(metadata == "Genotype")
  print(data)
  
  # load the annotation
  annotation <- read.delim(filepath_to_annotation,header=TRUE,row.names=1)
  annotation$feature <- row.names(annotation)
  annotation <- annotation %>% select(c("feature", "description"))
  annotation$feature <- gsub("-", ".", annotation$feature)
  print(annotation)
  # merge the data and annotation
  data <- merge(data, annotation, by = "feature")
  print(data)
  
  # group features into categories
  
  # subset the data based on the value column
  data_mut_significant <- data %>% filter(value == "MUT")
  data_het_significant <- data %>% filter(value == "HET")
  
  # return the subsetted data
  return(list(data_mut_significant, data_het_significant))
}

## Luminal comparisons --
tri_data <- load_and_process_data(filepath_to_significant_results = 
                                    "Trios/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv",
                                  filepath_to_annotation = 
                                    "Trios/differential_Pathway/annotated_pwy.tsv")

lt_data <- load_and_process_data(filepath_to_significant_results = 
                                   "Long_Term/differential_Pathway/PICRUST2_PWY_Luminal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv",
                                 filepath_to_annotation = 
                                   "Long_Term/differential_Pathway/annotated_pathway.tsv")

bl_data <- load_and_process_data(filepath_to_significant_results = 
                                     "Baseline/differential_Pathway/Baseline_JAX_PWY_Maaslin2_Sex_Genotype/significant_results.tsv",
                                   filepath_to_annotation = 
                                     "Baseline/differential_Pathway/annotated_pwy.tsv")

# combine the Luminal  data
trios <- as.data.frame(tri_data[[1]])
trios$Dataset <- c("Trios")
long_term <-as.data.frame(lt_data[[1]])
long_term$Dataset <- c("Long_Term")
baseline <- as.data.frame(bl_data[[1]])
baseline$Dataset <- c("Baseline_JAX")

all_ec_mut_significant <- rbind(trios, long_term,baseline)

# find intersecting features
v1 <- intersect(baseline$feature, long_term$feature) #6 features
v2 <- intersect(baseline$feature, trios$feature) #1 features
v3 <- intersect(trios$feature, long_term$feature) # 0 feature

v4 <- c(v1, v2,v3)
target <- unique(v4)

# subset the data based on the intersecting features
subset_all_ec_mut_significant <- all_ec_mut_significant[all_ec_mut_significant$feature %in% target, ]
subset_all_ec_het_significant <- all_ec_het_significant[all_ec_het_significant$feature %in% target, ]

df <- subset_all_ec_mut_significant

# Create the heatmap
ggplot(df, aes(x = Dataset, y = description, fill = coef)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  labs(x = "Dataset", y = "Feature") +
  theme_bw() +
  theme_cowplot(16)+
  ggtitle("Luminal Enriched MUT (yellow) vs Enriched WT (purple): 2 out of 3 datasets")

## Mucosal comparisons --
tri_data <- load_and_process_data(filepath_to_significant_results = 
                                    "Trios/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_SequencingRun_Site_Sex_Genotype_1-MouseID/significant_results.tsv",
                                  filepath_to_annotation = 
                                    "Trios/differential_Pathway/annotated_pwy.tsv")

lt_data <- load_and_process_data(filepath_to_significant_results = 
                                   "Long_Term/differential_Pathway/PICRUST2_PWY_Mucosal_Colon_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv",
                                 filepath_to_annotation = 
                                   "Long_Term/differential_Pathway/annotated_pathway.tsv")


# combine the Luminal  data
trios <- as.data.frame(tri_data[[1]])
trios$Dataset <- c("Trios")
long_term <-as.data.frame(lt_data[[1]])
long_term$Dataset <- c("Long_Term")

all_ec_mut_significant <- rbind(trios, long_term)

# find intersecting features
v3 <- intersect(trios$feature, long_term$feature) # 0 feature

