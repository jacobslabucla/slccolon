### Trios - Dataset wrangling ---
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(Maaslin2) #v 1.2.0
library(funrar)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_L6_Maaslin2.R")

metadata <- read.table("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(lumcol_meta) <- lumcol_meta$SampleID
lumcol <- lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(muccol_meta) <- muccol_meta$SampleID
muccol <- muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))


## Read in L6 aggregated counts -- 
lc_genus <- read.delim("Trios/collapsed_ASV/export_L",header=TRUE, row.names=1)
lc_genus <- lc_genus %>% select(names(c(lumcol_counts)))

mc_genus <- read.delim("Trios/collapsed_ASV/export_L6_Mucosal_Colon_notax_Trios_ASV_table_Silva_v138_1/feature-table.tsv",header=TRUE, row.names=1)
mc_genus <- mc_genus %>% select(names(c(muccol_counts)))

run_Maaslin2 <- function(counts_filepath, metadata_filepath, subset_string, fixed_effects_vector, results_string, random_effects_vector) {
  
  #input_data <- read.csv(counts_filepath, header=TRUE, row.names=1) # choose filtered non rarefied csv file
  
  df_input_data <- as.data.frame(counts_filepath)
  #df_input_data <- select(df_input_data, -c("taxonomy"))
  #df_input_data <- lumcol_counts
  
  transposed_input_data <- t(df_input_data)
  transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
  df_relative_ASV <- make_relative(transposed_input_data)
  df_relative_ASV <- as.data.frame(df_relative_ASV)
  Relative_Abundance <- summarize_all(df_relative_ASV, mean)
  Relative_Abundance <- as.data.frame(t(Relative_Abundance))
  
  readr::write_rds(Relative_Abundance,paste0("Trios/",subset_string,"Relative_Abundance-ASV.RDS"))
  
  #input_metadata <-read.delim(metadata_filepath,sep="\t",header=TRUE, row.names=1)
  input_metadata <- as.data.frame(metadata_filepath)
  #input_metadata <- lumcol_meta
  
  target <- colnames(df_input_data)
  input_metadata = input_metadata[match(target, row.names(input_metadata)),]
  target == row.names(input_metadata)
  
  df_input_metadata<-input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon", "Cecum"))
  df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels = c("Male", "Female"))
  df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
  df_input_metadata$Litter <- factor(df_input_metadata$Litter)
  df_input_metadata$SampleType <- factor(df_input_metadata$SampleType)
  df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
  sapply(df_input_metadata,levels)
  
  ?Maaslin2
  fit_data = Maaslin2(input_data=df_input_data, 
                      input_metadata=df_input_metadata, 
                      output = paste0("Trios/",subset_string, {{results_string}}), 
                      fixed_effects = {{fixed_effects_vector}},normalization="TSS",
                      random_effects= {{random_effects_vector}},
                      min_prevalence = 0.15,
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}

## Genus level Model 8 -- sus

# Luminal Colon-- Model fails to converge
fixed_effects <- c( "Site", "Sex", "Genotype")
random_effects <- c("Litter", "MouseID")
run_Maaslin2(lc_genus,lumcol_meta,"differential_taxa/L6_Luminal_Colon", 
             fixed_effects, "_L6_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter", random_effects)

# Mucosal Colon-- is ok 
fixed_effects <- c("Site", "Sex","Genotype")
random_effects <- c("Litter","MouseID")
run_Maaslin2(mc_genus,muccol_meta,"differential_taxa/L6_Mucosal_Colon", 
             fixed_effects, "_L6_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter", random_effects)

## Genus level Model 9 -- all good, no erros 

# Luminal Colon-- Model fails to converge
fixed_effects <- c("Sequencing_Run","Site", "Sex", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lc_genus,lumcol_meta,"differential_taxa/L6_Luminal_Colon", 
             fixed_effects, "_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID", random_effects)

# Mucosal Colon-- is ok 
fixed_effects <- c("Sequencing_Run","Site", "Sex", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2(mc_genus,muccol_meta,"differential_taxa/L6_Mucosal_Colon", 
             fixed_effects, "_L6_Maaslin2_Sequencing_Run_Site_Sex_Genotype_1-MouseID", random_effects)
