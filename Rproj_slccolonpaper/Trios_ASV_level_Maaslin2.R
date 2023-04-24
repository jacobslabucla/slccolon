### Trios - Dataset wrangling ---
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(Maaslin2) #v 1.2.0
library(funrar)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_ASV_level_Maaslin2.R")

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

# Cecum
cecum_meta <- metadata %>% filter(Site=="Cecum", SampleID %in% names(counts))
row.names(cecum_meta) <- cecum_meta$SampleID
cecum <- cecum_meta$SampleID
cecum_counts <- counts %>% select(all_of(cecum))

# Proximal_Colon
pc_meta <- metadata %>% filter(Site=="Proximal_Colon", SampleID %in% names(counts))
row.names(pc_meta) <- pc_meta$SampleID
pc <- pc_meta$SampleID
pc_counts <- counts %>% select(all_of(pc))

# Distal_Colon
dc_meta <- metadata %>% filter(Site=="Distal_Colon", SampleID %in% names(counts))
row.names(dc_meta) <- dc_meta$SampleID
dc <- dc_meta$SampleID
dc_counts <- counts %>% select(all_of(dc))

# Cecum Luminal 
lum_cecum_meta <- metadata %>% filter(Site=="Cecum"& SampleType=="Luminal", SampleID %in% names(counts))
row.names(lum_cecum_meta) <- lum_cecum_meta$SampleID
lum_cecum <- lum_cecum_meta$SampleID
lum_cecum_counts <- counts %>% select(all_of(lum_cecum))

# Proximal_Colon Luminal
lum_pc_meta <- metadata %>% filter(Site=="Proximal_Colon"& SampleType=="Luminal", SampleID %in% names(counts))
row.names(lum_pc_meta) <- lum_pc_meta$SampleID
lum_pc <- lum_pc_meta$SampleID
lum_pc_counts <- counts %>% select(all_of(lum_pc))

# Distal_Colon Luminal 
lum_dc_meta <- metadata %>% filter(Site=="Distal_Colon"& SampleType=="Luminal", SampleID %in% names(counts))
row.names(lum_dc_meta) <- lum_dc_meta$SampleID
lum_dc <- lum_dc_meta$SampleID
lum_dc_counts <- counts %>% select(all_of(lum_dc))


## Run Maaslin2 with prevalence filtering to 15% on subsets -- 

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
  
  readr::write_rds(Relative_Abundance,paste0("Trios/Relative_Abundance-",subset_string,"-ASV.RDS"))
  
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

## Model 1

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Sex", "Site","Genotype")
random_effects <- c("Litter", "MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-Litter_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Sex", "Site","Genotype")
random_effects <- c("Litter", "MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-Litter_1-MouseID", random_effects)

## Model 2

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Sex", "Litter","Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Litter_Site_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Sex", "Litter","Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Litter_Site_Genotype_1-MouseID", random_effects)

## Model 3

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Sex", "Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Sex","Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID", random_effects)

## Model 4 - no errors 

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Litter", "Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Site_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Litter","Site","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Site_Genotype_1-MouseID", random_effects)

## Model 5 - no errors 

# Luminal Colon
fixed_effects <- c("Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_Genotype_1-MouseID", random_effects)

## Model 6 - rank deficient 

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Litter", "Sex","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Sex_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Litter","Sex","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Sex_Genotype_1-MouseID", random_effects)

## Model 7 - rank deficient 

# Luminal Colon
fixed_effects <- c("Sequencing_Run", "Litter","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Genotype_1-MouseID", random_effects)

# Mucosal Colon
fixed_effects <- c("Sequencing_Run", "Litter","Genotype")
random_effects <- c("MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Litter_Genotype_1-MouseID", random_effects)

## Model 8 -  no errors

# Luminal Colon
fixed_effects <- c( "Site", "Sex", "Genotype")
random_effects <- c("Litter", "MouseID")
run_Maaslin2(lumcol_counts,lumcol_meta,"Luminal Colon", fixed_effects, "_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter", random_effects)

# Mucosal Colon
fixed_effects <- c("Site", "Sex","Genotype")
random_effects <- c("Litter","MouseID")
run_Maaslin2(muccol_counts,muccol_meta,"Mucosal Colon", fixed_effects, "_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter", random_effects)

## Run Maaslin2 with prevalence filtering to 15% on each site dataset -- 

# Cecum
fixed_effects <- c("Sequencing_Run", "Sex", "Litter", "Genotype")
random_effects <- c("")
run_Maaslin2(cecum_counts,cecum_meta,"Cecum", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Litter_Genotype", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "Genotype")
random_effects <- c("Litter")
run_Maaslin2(cecum_counts,cecum_meta,"Cecum", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "SampleType","Genotype")
random_effects <- c("Litter")
run_Maaslin2(cecum_counts,cecum_meta,"Cecum", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex","Genotype")
random_effects <- c("Litter")
run_Maaslin2(lum_cecum_counts,lum_cecum_meta,"Lum_Cecum", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter", random_effects)


# Proximal Colon 
fixed_effects <- c("Sequencing_Run", "Sex", "Litter", "Genotype")
random_effects <- c("")
run_Maaslin2(pc_counts,pc_meta,"Proximal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Litter_Genotype", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "Genotype")
random_effects <- c("Litter")
run_Maaslin2(pc_counts,pc_meta,"Proximal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "SampleType","Genotype")
random_effects <- c("Litter")
run_Maaslin2(pc_counts,pc_meta,"Proximal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "Genotype")
random_effects <- c("Litter")
run_Maaslin2(lum_pc_counts,lum_pc_meta,"Lum_Proximal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter", random_effects)

# Distal Colon 
fixed_effects <- c("Sequencing_Run", "Sex", "Litter", "Genotype")
random_effects <- c("")
run_Maaslin2(dc_counts,dc_meta,"Distal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Litter_Genotype", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "Genotype")
random_effects <- c("Litter")
run_Maaslin2(dc_counts,dc_meta,"Distal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "SampleType","Genotype")
random_effects <- c("Litter")
run_Maaslin2(dc_counts,dc_meta,"Distal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter", random_effects)

fixed_effects <- c("Sequencing_Run", "Sex", "Genotype")
random_effects <- c("Litter")
run_Maaslin2(lum_dc_counts,lum_dc_meta,"Lum_Distal_Colon", fixed_effects, "_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter", random_effects)


## Add annotations to Maaslin2 results --

## Subset Results for Model 3 

# Luminal Colon: 0 significant features by Genotype 

data<-read.table("Trios/Luminal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Luminal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/annotated_significant_results.csv")

# Mucosal Colon: 1 significant features by Genotype 

data<-read.table("Trios/Mucosal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Mucosal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/annotated_significant_results.csv")

## Subset Results for Model 5 

# Luminal Colon: 0 significant features by Genotype 

data<-read.table("Trios/Luminal Colon_ASV_Maaslin2_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
#write.csv(data, "Trios/Luminal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/annotated_significant_results.csv")

# Mucosal Colon: 0 significant features by Genotype 

data<-read.table("Trios/Mucosal Colon_ASV_Maaslin2_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
#write.csv(data, "Trios/Mucosal Colon_ASV_Maaslin2_SeqRun_Sex_Site_Genotype_1-MouseID/annotated_significant_results.csv")

## Subset Results for Model 8 

# Luminal Colon: 3 significant features by Genotype at q<0.05

data<-read.table("Trios/Luminal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Luminal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/annotated_significant_results.csv")

# Mucosal Colon: 1 significant features by Genotype at q<0.05

data<-read.table("Trios/Mucosal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/significant_results.tsv", header=TRUE)
data <- data  %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Mucosal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/annotated_significant_results.csv") 

# Proximal Colon 

data<-read.table("Trios/Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Cecum 
data<-read.table("Trios/Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Divided by Site -
# Distal Colon 

data<-read.table("Trios/Distal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Distal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Proximal Colon 

data<-read.table("Trios/Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Cecum 
data<-read.table("Trios/Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05) %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Lum Cecum - 5 features with q <0.25
data<-read.table("Trios/Lum_Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Lum_Cecum_ASV_Maaslin2_SeqRun_Sex_Type_Genotype_1-Litter/annotated_significant_results.csv")

# Lum Proximal Colon - 2 features with q <0.25
data<-read.table("Trios/Lum_Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>%  filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Lum_Proximal_Colon_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter/annotated_significant_results.csv")

# Lum Distal Colon - 2 features with q <0.25
data<-read.table("Trios/Lum_Distal_Colon_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(metadata=="Genotype")

data<- (merge(data, annotation, by = 'feature'))
write.csv(data, "Trios/Lum_Distal_Colon_ASV_Maaslin2_SeqRun_Sex_Genotype_1-Litter/annotated_significant_results.csv")


