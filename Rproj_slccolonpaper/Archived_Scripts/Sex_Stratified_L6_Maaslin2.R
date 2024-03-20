### Trios - Dataset wrangling ---
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(Maaslin2) #v 1.2.0
library(funrar)
library(plyr)

here::i_am("Rproj_slccolonpaper/Archived_Scripts/Sex_Stratified_L6_Maaslin2.R")

metadata <- read.table(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"), header=TRUE)
counts <- read.table(here("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv"), header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon M
m_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts),Sex=="Male")
row.names(m_lumcol_meta) <- m_lumcol_meta$SampleID
m_lumcol <- m_lumcol_meta$SampleID
m_lumcol_counts <- counts %>% select(all_of(m_lumcol))

# Luminal Colon F
f_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts),Sex=="Female")
row.names(f_lumcol_meta) <- f_lumcol_meta$SampleID
f_lumcol <- f_lumcol_meta$SampleID
f_lumcol_counts <- counts %>% select(all_of(f_lumcol))

# Mucosal Colon F
f_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts), Sex=="Female")
row.names(f_muccol_meta) <- f_muccol_meta$SampleID
f_muccol <- f_muccol_meta$SampleID
f_muccol_counts <- counts %>% select(all_of(f_muccol))

# Mucosal Colon M
m_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts), Sex=="Male")
row.names(m_muccol_meta) <- m_muccol_meta$SampleID
m_muccol <- m_muccol_meta$SampleID
m_muccol_counts <- counts %>% select(all_of(m_muccol))

## Read in L6 aggregated counts -- 
lc_genus <- read.delim("Trios/collapsed_ASV/export_L6_Luminal_Colon_notax_Trios_ASV_table_Silva_v138_1/feature-table.tsv",header=TRUE, row.names=1)
f_lc_genus <- lc_genus %>% select(names(c(f_lumcol_counts)))
m_lc_genus <- lc_genus %>% select(names(c(m_lumcol_counts)))

mc_genus <- read.delim("Trios/collapsed_ASV/export_L6_Mucosal_Colon_notax_Trios_ASV_table_Silva_v138_1/feature-table.tsv",header=TRUE, row.names=1)
f_mc_genus <- mc_genus %>% select(names(c(f_muccol_counts)))
m_mc_genus <- mc_genus %>% select(names(c(m_muccol_counts)))

run_Maaslin2_stratified <- function(counts_filepath, metadata_filepath, subset_string, fixed_effects_vector, results_string, random_effects_vector,
                                    ref_string) {
  
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
  print(target == row.names(input_metadata))
  
  df_input_metadata<-input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon", "Cecum"))
  df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels = c("Male", "Female"))
  df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
  df_input_metadata$Litter <- factor(df_input_metadata$Litter)
  df_input_metadata$SampleType <- factor(df_input_metadata$SampleType)
  df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, levels=c("WT","HET", "MUT"))
  #df_input_metadata$Genotype <- revalue( df_input_metadata$Genotype, c("WT"="a_WT", "HET" = "b_HET", "MUT" ="c_MUT"))
  sapply(df_input_metadata,levels)
  
  ?Maaslin2
  fit_data = Maaslin2(input_data=df_input_data, 
                      input_metadata=df_input_metadata, 
                      output = paste0("Trios/",subset_string, {{results_string}}), 
                      fixed_effects = {{fixed_effects_vector}},normalization="TSS",
                      random_effects= {{random_effects_vector}},
                      min_prevalence = 0.15,
                      reference={{ref_string}},
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}

# Luminal Colon-- neither have convergence errors
fixed_effects <- c("Litter","Site", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2_stratified(f_lc_genus,f_lumcol_meta,"differential_taxa/F_L6_Luminal_Colon", 
             fixed_effects, "_L6_Maaslin2_Litter_Site_Genotype_1-MouseID", random_effects,
             c(("Litter,one"), ("Site,Distal_Colon"),("Genotype,WT")))
run_Maaslin2_stratified(m_lc_genus,m_lumcol_meta,"differential_taxa/M_L6_Luminal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Litter_Site_Genotype_1-MouseID", random_effects,
                        c("Litter,one", "Site,Distal_Colon","Genotype,WT"))

fixed_effects <- c("Site", "Genotype")
run_Maaslin2_stratified(f_lc_genus,f_lumcol_meta,"differential_taxa/F_L6_Luminal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Site_Genotype_1-MouseID", random_effects,
                        c("Site,Distal_Colon","Genotype,WT"))
run_Maaslin2_stratified(m_lc_genus,m_lumcol_meta,"differential_taxa/M_L6_Luminal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Site_Genotype_1-MouseID", random_effects,
                        c("Site,Distal_Colon","Genotype,WT"))

fixed_effects <- c("Sequencing_Run","Site", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2_stratified(m_lc_genus,m_lumcol_meta,"differential_taxa/M_L6_Luminal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Sequencing_Run_Site_Genotype_1-MouseID", random_effects,
                        c("Sequencing_Run,SLC","Site,Distal_Colon","Genotype,WT"))

# Mucosal Colon-- is ok 
fixed_effects <- c("Litter","Site", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2_stratified(f_mc_genus,f_muccol_meta,"differential_taxa/F_L6_Mucosal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Litter_Site_Genotype_1-MouseID", random_effects,
                        c("Litter,one", "Site,Distal_Colon","Genotype,WT"))
run_Maaslin2_stratified(m_mc_genus,m_muccol_meta,"differential_taxa/M_L6_Mucosal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Litter_Site_Genotype_1-MouseID", random_effects,
                        c("Litter,one", "Site,Distal_Colon","Genotype,WT"))

fixed_effects <- c("Site", "Genotype")
run_Maaslin2_stratified(f_mc_genus,f_muccol_meta,"differential_taxa/F_L6_Mucosal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Site_Genotype_1-MouseID", random_effects,
                        c("Site,Distal_Colon","Genotype,WT"))
run_Maaslin2_stratified(m_mc_genus,m_muccol_meta,"differential_taxa/M_L6_Mucosal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Site_Genotype_1-MouseID", random_effects,
                        c("Site,Distal_Colon","Genotype,WT"))

fixed_effects <- c("Sequencing_Run","Site", "Genotype")
random_effects <- c("MouseID")
run_Maaslin2_stratified(m_mc_genus,m_muccol_meta,"differential_taxa/M_L6_Mucosal_Colon", 
                        fixed_effects, "_L6_Maaslin2_Sequencing_Run_Site_Genotype_1-MouseID", random_effects,
                        c("Sequencing_Run,SLC","Site,Distal_Colon","Genotype,WT"))

### Fecal Pellet data ---
metadata <- read.table(here("Baseline/starting_files/Baseline_Metadata.tsv"), header=TRUE)
counts <- read.table(here("Baseline/starting_files/Baseline_ASV_table_Silva_v138_1.tsv"), header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Get counts only from the JAX mice

f_JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts), Sex=="Female")
row.names(f_JAX_meta) <- f_JAX_meta$SampleID
f_JAX <- f_JAX_meta$SampleID
f_JAX_counts <- counts %>% select(all_of(f_JAX))

m_JAX_meta <- metadata %>% filter(Background=="JAX", SampleID %in% names(counts), Sex=="Male")
row.names(m_JAX_meta) <- m_JAX_meta$SampleID
m_JAX <- m_JAX_meta$SampleID
m_JAX_counts <- counts %>% select(all_of(m_JAX))

## Read in L6 aggregated counts -- 
lc_genus <- read.delim("Baseline/collapsed_ASV/export_L6_JAX_Baseline_ASV_table_Silva_v138_1/feature-table.tsv",header=TRUE, row.names=1)
f_lc_genus <- lc_genus %>% select(names(c(f_JAX_counts)))
m_lc_genus <- lc_genus %>% select(names(c(m_JAX_counts)))

### Run Maaslin2 and get table of relative abundances

run_Maaslin2_fp <- function(counts, metadata, subset_string,ref_string) {
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
                      output = paste0("Baseline/",subset_string,"_Maaslin2_Genotype"), 
                      fixed_effects = c("Genotype"),normalization="TSS", 
                      min_prevalence = 0.15,
                      reference={{ref_string}},
                      transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)
}

run_Maaslin2_fp(f_lc_genus,f_JAX_meta, "differential_taxa/f_JAX_L6", "Genotype,WT")

run_Maaslin2_fp(m_lc_genus,m_JAX_meta, "differential_taxa/m_JAX_L6", "Genotype,WT")
