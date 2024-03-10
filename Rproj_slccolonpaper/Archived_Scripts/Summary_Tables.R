library(pivottabler)
library(dplyr)
library(knitr)

## Feed in Metadata for the 3 cohorts --

metadata <- read.csv(here("Long_Term/starting_files/SLC_LT_metadata.csv"), header=TRUE)

# Make a pivot table 
metadata$MouseID <- factor(metadata$MouseID)
slt_summary <- metadata %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

slt_summary %>% kable

## Feed in Metadata for the 3 cohorts --

metadata <- read.delim(here("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv"), header=TRUE)

# Make a pivot table 
metadata$MouseID <- factor(metadata$MouseID)
trios_summary <- metadata %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

trios_summary %>% kable


# Feed in Metadata for the 3 cohorts --

metadata <- read.csv("Baseline/starting_files/Baseline_Metadata.csv", header=TRUE)
sapply(metadata,levels)

# Make a pivot table 
metadata$MouseID <- factor(metadata$MouseID)
baseline_summary <- metadata %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

baseline_summary %>% kable

pff_mptp <- metadata %>% filter(Line!="Q22")

# Make a pivot table 

baseline_summary <- pff_mptp %>%
  group_by(Sex, Genotype) %>%
  summarize(MouseID = n_distinct(MouseID)) 

baseline_summary %>% kable
