library(vegan) #v 2.6.2

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files
here::i_am("Rproj_slccolonpaper/Trios_RS_Jensen_Shannon.R")

metadata <- read.table("Trios/starting_files/SLC_TOTAL_OCT2020_FULL_Metadata.tsv", header=TRUE)
counts <- read.table("Trios/starting_files/Trios_ASV_table_Silva_v138_1.tsv", header = TRUE, row.names=1)

## Store taxonomy in an annotation file --
annotation <- tibble::rownames_to_column(counts, "feature") %>% select(c("feature", "taxonomy"))
counts <- counts %>% select(-c("taxonomy"))

## Apply minimum sequencing depth threshold --
counts <- counts[colSums(counts) >= 10000]

## Split counts into colon subsets -- 

# Luminal Colon 
trios_lumcol_meta <- metadata %>% filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(trios_lumcol_meta) <- trios_lumcol_meta$SampleID
lumcol <- trios_lumcol_meta$SampleID
lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon
trios_muccol_meta <- metadata %>% filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(trios_muccol_meta) <- trios_muccol_meta$SampleID
muccol <- trios_muccol_meta$SampleID
muccol_counts <- counts %>% select(all_of(muccol))

# Luminal Colon no HET 
nohet_trios_lumcol_meta <- metadata %>%
  filter(Genotype!="HET")%>%
  filter(Subset=="Luminal_Colon", SampleID %in% names(counts))
row.names(nohet_trios_lumcol_meta) <- nohet_trios_lumcol_meta$SampleID
lumcol <- nohet_trios_lumcol_meta$SampleID
nohet_lumcol_counts <- counts %>% select(all_of(lumcol))

# Mucosal Colon no HET
nohet_trios_muccol_meta <- metadata %>% 
  filter(Genotype!="HET")%>%
  filter(Subset=="Mucosal_Colon", SampleID %in% names(counts))
row.names(nohet_trios_muccol_meta) <- nohet_trios_muccol_meta$SampleID
muccol <- nohet_trios_muccol_meta$SampleID
nohet_muccol_counts <- counts %>% select(all_of(muccol))

## Prevalence filter datasets -- 
# Luminal Colon
trios_lumcol_counts <- prevalence_filter(lumcol_counts,13)

# Mucosal Colon 
trios_muccol_counts <- prevalence_filter(muccol_counts,12)

# Luminal Colon - no HET 
0.15*60 #9 samples
nohet_trios_lumcol_counts <- prevalence_filter(nohet_lumcol_counts,9)

# Mucosal Colon - no HET 
0.15*52 #8 samples
nohet_trios_muccol_counts <- prevalence_filter(nohet_muccol_counts,8)

## Calculate RS Jensen Shannon distance matrix -- 


trios_muccol.dist <- calculate_rsjensen(trios_muccol_counts)
trios_lumcol.dist <- calculate_rsjensen(trios_lumcol_counts)
nohet_trios_muccol.dist <- calculate_rsjensen(nohet_trios_muccol_counts)
nohet_trios_lumcol.dist <- calculate_rsjensen(nohet_trios_lumcol_counts)




## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=lumcol.dist,
                    counts = trios_lumcol_counts,
                    metadata = trios_lumcol_meta,
                    title="Trios - Luminal Colon RS Jensen",
                    colorvariable = Genotype,
                    colorvector = cols,
                    wa_scores_filepath = "Trios/beta_diversity/LumCol_Top_Taxa_PcoA.csv")

trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=muccol.dist,
                                     counts = trios_muccol_counts,
                                     metadata = trios_muccol_meta,
                                     title="Trios - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/MucCol_Top_Taxa_PcoA.csv")

nohet_trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=nohet_trios_lumcol.dist,
                                     counts = nohet_trios_lumcol_counts,
                                     metadata = nohet_trios_lumcol_meta,
                                     title="Trios - Luminal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/nohet_LumCol_Top_Taxa_PcoA.csv")

nohet_trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=nohet_trios_muccol.dist,
                                     counts = nohet_trios_muccol_counts,
                                     metadata = nohet_trios_muccol_meta,
                                     title="Trios - Mucosal Colon RS Jensen",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/nohet_MucCol_Top_Taxa_PcoA.csv")

plot_grid(trios_lc_pcoa, trios_mc_pcoa,labels=c("A","B"),label_size=20)
plot_grid(nohet_trios_lc_pcoa, nohet_trios_mc_pcoa,labels=c("A","B"),label_size=20)

## Statistics --

## PERMANOVA

# Luminal Colon
data.dist<-trios_lumcol.dist
metadata <- trios_lumcol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon
data.dist<-trios_muccol.dist
metadata <- trios_muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Luminal Colon - no HETs
data.dist<-nohet_trios_lumcol.dist
metadata <- nohet_trios_lumcol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

# Mucosal Colon - no HETs
data.dist<-nohet_trios_muccol.dist
metadata <- nohet_trios_muccol_meta

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Litter + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Genotype, data=metadata, permutations=10000)
data.adonis$aov.tab



## Repeat-Measures-Aware 
# Luminal Colon
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_lumcol.dist,
                       trios_lumcol_meta,
                       site,
                       mouseID,
                       order_vector)

site <- c("Site")
mouseID <- c("Sequencing_Run","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_lumcol.dist,
                       trios_lumcol_meta,
                       site,
                       mouseID,
                       order_vector)

names(nohet_trios_lumcol_meta)
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(nohet_trios_lumcol.dist,
                       nohet_trios_lumcol_meta,
                       site,
                       mouseID,
                       order_vector)

# Mucosal Colon
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_muccol.dist,
                       trios_muccol_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sequencing_Run","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_muccol.dist,
                       trios_muccol_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(nohet_trios_muccol.dist,
                       nohet_trios_muccol_meta,
                       site,
                       mouseID,
                       order_vector)

