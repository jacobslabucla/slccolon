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

# Luminal SI
trios_lumsi_meta <- metadata %>% filter(Subset=="Luminal_SI", SampleID %in% names(counts))
row.names(trios_lumsi_meta) <- trios_lumsi_meta$SampleID
lumsi <- trios_lumsi_meta$SampleID
lumsi_counts <- counts %>% select(all_of(lumsi))

# Mucosal sion
trios_mucsi_meta <- metadata %>% filter(Subset=="Mucosal_SI", SampleID %in% names(counts))
row.names(trios_mucsi_meta) <- trios_mucsi_meta$SampleID
mucsi <- trios_mucsi_meta$SampleID
mucsi_counts <- counts %>% select(all_of(mucsi))

# Luminal sion no HET 
nohet_trios_lumsi_meta <- metadata %>%
  filter(Genotype!="HET")%>%
  filter(Subset=="Luminal_SI", SampleID %in% names(counts))
row.names(nohet_trios_lumsi_meta) <- nohet_trios_lumsi_meta$SampleID
lumsi <- nohet_trios_lumsi_meta$SampleID
nohet_lumsi_counts <- counts %>% select(all_of(lumsi))

# Mucosal sion no HET
nohet_trios_mucsi_meta <- metadata %>% 
  filter(Genotype!="HET")%>%
  filter(Subset=="Mucosal_SI", SampleID %in% names(counts))
row.names(nohet_trios_mucsi_meta) <- nohet_trios_mucsi_meta$SampleID
mucsi <- nohet_trios_mucsi_meta$SampleID
nohet_mucsi_counts <- counts %>% select(all_of(mucsi))

## Prevalence filter datasets -- 
# Luminal SI
88*0.15
trios_lumsi_counts <- prevalence_filter(lumsi_counts,13)

# Mucosal sion 
81*0.15
trios_mucsi_counts <- prevalence_filter(mucsi_counts,12)

# Luminal SI - no HET 
0.15*60 #9 samples
nohet_trios_lumsi_counts <- prevalence_filter(nohet_lumsi_counts,9)

# Mucosal SI- no HET 
0.15*52 #8 samples
nohet_trios_mucsi_counts <- prevalence_filter(nohet_mucsi_counts,8)

## Calculate RS Jensen Shannon distance matrix -- 


trios_mucsi.dist <- calculate_rsjensen(trios_mucsi_counts)
trios_lumsi.dist <- calculate_rsjensen(trios_lumsi_counts)
nohet_trios_mucsi.dist <- calculate_rsjensen(nohet_trios_mucsi_counts)
nohet_trios_lumsi.dist <- calculate_rsjensen(nohet_trios_lumsi_counts)




## Principal Coordinates Analysis -- 

cols <- c("WT"="black", "HET"= "blue", "MUT"="red")

trios_lsi_pcoa <- generate_pcoA_plots(distance_matrix=trios_lumsi.dist,
                                     counts = trios_lumsi_counts,
                                     metadata = trios_lumsi_meta,
                                     title="Luminal SI",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/LumSI_Top_Taxa_PcoA.csv")

trios_msi_pcoa <- generate_pcoA_plots(distance_matrix=trios_mucsi.dist,
                                     counts = trios_mucsi_counts,
                                     metadata = trios_mucsi_meta,
                                     title="Mucosal SI",
                                     colorvariable = Genotype,
                                     colorvector = cols,
                                     wa_scores_filepath = "Trios/beta_diversity/mucsi_Top_Taxa_PcoA.csv")

nohet_trios_lc_pcoa <- generate_pcoA_plots(distance_matrix=nohet_trios_lumsi.dist,
                                           counts = nohet_trios_lumsi_counts,
                                           metadata = nohet_trios_lumsi_meta,
                                           title="Trios - Luminal Colon RS Jensen",
                                           colorvariable = Genotype,
                                           colorvector = cols,
                                           wa_scores_filepath = "Trios/beta_diversity/nohet_lumsi_Top_Taxa_PcoA.csv")

nohet_trios_mc_pcoa <- generate_pcoA_plots(distance_matrix=nohet_trios_mucsi.dist,
                                           counts = nohet_trios_mucsi_counts,
                                           metadata = nohet_trios_mucsi_meta,
                                           title="Trios - Mucosal Colon RS Jensen",
                                           colorvariable = Genotype,
                                           colorvector = cols,
                                           wa_scores_filepath = "Trios/beta_diversity/nohet_mucsi_Top_Taxa_PcoA.csv")

plot_grid(trios_lsi_pcoa, trios_msi_pcoa,labels=c("A","B"),label_size=20)
plot_grid(nohet_trios_lc_pcoa, nohet_trios_mc_pcoa,labels=c("A","B"),label_size=20)

## Statistics --

## PERMANOVA

# Luminal SI
data.dist<-trios_lumsi.dist
metadata <- trios_lumsi_meta

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
data.dist<-trios_mucsi.dist
metadata <- trios_mucsi_meta

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
data.dist<-nohet_trios_lumsi.dist
metadata <- nohet_trios_lumsi_meta

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
data.dist<-nohet_trios_mucsi.dist
metadata <- nohet_trios_mucsi_meta

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
run_repeated_PERMANOVA(trios_lumsi.dist,
                       trios_lumsi_meta,
                       site,
                       mouseID,
                       order_vector)

site <- c("Site")
mouseID <- c("Sequencing_Run","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_lumsi.dist,
                       trios_lumsi_meta,
                       site,
                       mouseID,
                       order_vector)

names(nohet_trios_lumsi_meta)
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(nohet_trios_lumsi.dist,
                       nohet_trios_lumsi_meta,
                       site,
                       mouseID,
                       order_vector)

# Mucosal Colon
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_mucsi.dist,
                       trios_mucsi_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sequencing_Run","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Sex","Site","Genotype")
run_repeated_PERMANOVA(trios_mucsi.dist,
                       trios_mucsi_meta,
                       site,
                       mouseID,
                       order_vector)
site <- c("Site")
mouseID <- c("Sequencing_Run","Litter","Sex","Genotype","MouseID")
order_vector <- c("Sequencing_Run","Litter","Sex","Site","Genotype")
run_repeated_PERMANOVA(nohet_trios_mucsi.dist,
                       nohet_trios_mucsi_meta,
                       site,
                       mouseID,
                       order_vector)

