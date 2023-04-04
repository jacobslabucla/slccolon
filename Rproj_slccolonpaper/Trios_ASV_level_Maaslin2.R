### Trios - Dataset wrangling ---
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.1.0
library(Maaslin2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Maaslin2")

setwd("/home/julianne/Documents/slccolonpaper/slccolon/")
here::i_am("Rproj_slccolonpaper/Trios_ASV_level_Maaslin2.R")

metadata <- read.table("Trios/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "ASV") %>% select(c("ASV", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
lumcol <- lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
muccol <- muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))

## Run Maaslin2 -- 
