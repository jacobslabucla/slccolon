### Trios - Taxa summary plots 
### Date : 4/4/23
library(here) #v 1.0.1
library(dplyr) #v 1.0.7
library(Maaslin2) #v 1.2.0
library(funrar)

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_ASV_level_Maaslin2.R")

metadata <- read.table("Trios/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

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

## Convert counts into relative abundances --

#Taxa are columns, samples are rows
lumcol_counts<- as.matrix(t(lumcol_counts))
lumcol_counts<-funrar::make_relative(lumcol_counts)
sum(lumcol_counts[1,])

#Samples are columns, taxa are rows. Calculate average abundance of each taxa across all samples
lumcol_counts<-as.data.frame(t(lumcol_counts))
toptaxa<- rowMeans(lumcol_counts)
length(names(lumcol_counts))

lumcol_counts$averageRA <-toptaxa
lumcol_counts <- lumcol_counts %>%
  dplyr::mutate(keeptaxa = ifelse(averageRA >0.01, row.names(lumcol_counts), "Other"))
lumcol_counts <-select(lumcol_counts,-averageRA)

taxa<-lumcol_counts$keeptaxa
lumcol_counts <- select(lumcol_counts,-keeptaxa)
lumcol_counts <- as.matrix(sapply(lumcol_counts,as.numeric))
lumcol_counts <- as.data.frame(prop.table(lumcol_counts,2))

L2_lum$Taxa <-taxa
L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
L2_lum$Value <- L2_lum$Value * 100
