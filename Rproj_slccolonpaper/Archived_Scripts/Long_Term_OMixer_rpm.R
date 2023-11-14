library(omixerRpm)
library(here)
library(dplyr)

setwd("/home/julianne/Documents/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Long_Term_OMixer_rpm.R")

metadata <- read.table("Long_Term/starting_files/SLC_LT_metadata.tsv", header=TRUE)
counts <- read.table("Long_Term/starting_files/SLT_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)
pathway <- read.delim("Long_Term/starting_files/picrust2_output_SLT_ASV_table_Silva_v138_1.qza/export_ko_metagenome/feature-table.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))
pathway <- pathway %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Grab samples in the pathway table that are present in counts --
names(pathway)
pathway <- pathway %>% select(c(names(counts)))

#load the database (select between gut-brain and gut-metabolic)
listDB()

# Gut Metabolic
db_metabolic<-loadDefaultDB() # by default is metabolic

#Gut Brain
db_brain <- loadDB("GBMs.v1.0")

#read table in special format 
dat <- pathway
dat <- tibble::rownames_to_column(dat, "KO")

#After utilizing GOMIXER to find optimal module coverage value
mods <- rpm(dat, minimum.coverage=0.6, annotation = 1, module.db=db_brain) #coverage- what percent of KO need to be present in one module for one mod to be present
mods <- rpm(dat, minimum.coverage=0.6, annotation = 1, module.db=db_metabolic) #coverage- what percent of KO need to be present in one module for one mod to be present


