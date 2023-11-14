library(Maaslin2)
library(funrar)
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Baseline_JAX_ASV_L6_Maaslin2.R")

metadata <- read.table("Baseline/starting_files/Baseline_Metadata.tsv", header=TRUE)
counts <- read.table("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Get counts only from the JAX mice

JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts))
row.names(JAX_meta) <- JAX_meta$SampleID
JAX <- JAX_meta$SampleID
JAX_counts <- counts %>% select(all_of(JAX))

## Read in L6 aggregated counts -- 
lc_genus <- read.delim("Baseline/collapsed_ASV/export_L6_JAX_Baseline_ASV_table_Silva_v138_1/feature-table.tsv",header=TRUE, row.names=1)
lc_genus <- lc_genus %>% select(names(c(JAX_counts)))

### Run Maaslin2 and get table of relative abundances

run_Maaslin2 <- function(counts, metadata, subset_string) {
  input_data <- as.data.frame(counts)
  df_input_data <- as.data.frame(input_data)
  
  transposed_input_data <- t(df_input_data)
  transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
  df_relative_ASV <- make_relative(transposed_input_data)
  df_relative_ASV <- as.data.frame(df_relative_ASV)
  Relative_Abundance <- summarize_all(df_relative_ASV, mean)
  Relative_Abundance <- as.data.frame(t(Relative_Abundance))
  
  readr::write_rds(Relative_Abundance,paste0("Baseline/",subset_string,"Relative_Abundance_ASV.RDS"))
  
  input_metadata <- as.data.frame(metadata)
  
  target <- colnames(df_input_data)
  input_metadata = input_metadata[match(target, row.names(input_metadata)),]
  target == row.names(input_metadata)
  
  df_input_metadata<-input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
  df_input_metadata$Sex <- factor(df_input_metadata$Sex)
  sapply(df_input_metadata,levels)
  
  ?Maaslin2
  fit_data = Maaslin2(input_data=df_input_data, 
                      input_metadata=df_input_metadata, 
                      output = paste0("Baseline/",subset_string,"_Maaslin2_Sex_Genotype"), 
                      fixed_effects = c("Sex","Genotype"),normalization="TSS", 
                      min_prevalence = 0.15,
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}

run_Maaslin2(JAX_counts,JAX_meta, "differential_taxa/JAX_ASV")

run_Maaslin2(lc_genus,JAX_meta, "differential_taxa/JAX_L6")

### Visualize Results ---

## ASV level : Make a taxa dotplot for SLC Trios q<0.25 --
data<-read.table("Baseline/differential_taxa/JAX_ASV_Maaslin2_Sex_Genotype/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

jax_mut <- data %>% filter(value=="MUT")
jax_het <- data %>% filter(value=="HET")

# Luminal Colon 
make_taxa_dotplot(jax_mut,
                  "Baseline/starting_files/Taxonomy_Key.tsv",
                  "Baseline/differential_taxa/JAX_ASVRelative_Abundance_ASV.RDS",
                  "Baseline JAX MUT: ASV~ Sex + Genotype q,0.25",
                  colorvariable = Phylum)


make_taxa_dotplot(jax_het,
                  "Baseline/starting_files/Taxonomy_Key.tsv",
                  "Baseline/differential_taxa/JAX_ASVRelative_Abundance_ASV.RDS",
                  "Baseline JAX HET: ASV~ Sex + Genotype q,0.25",
                  colorvariable = Phylum)

