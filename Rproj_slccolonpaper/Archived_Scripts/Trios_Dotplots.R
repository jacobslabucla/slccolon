
library(dada2)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files


## Make a taxa dotplot for SLC Trios q<0.25 --
data<-read.table("Trios/Luminal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")

# Luminal Colon 
make_taxa_dotplot(data,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "Trios Luminal Colon: ASV ~ Site + Sex + Genotype + (1|MouseID) + (1|Litter) q<0.25")

make_taxa_dotplot(trios_mut,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "Trios Luminal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) + (1|Litter) q<0.25")

make_taxa_dotplot(trios_het,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "Trios Luminal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) + (1|Litter) q<0.25")

# Mucosal Colon
data<-read.table("Trios/Mucosal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID_1-Litter/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")

make_taxa_dotplot(data,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "Trios Mucosal Colon: ASV ~ Site + Sex + Genotype +(1|MouseID) + (1|Litter) q<0.25")

make_taxa_dotplot(trios_mut,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "Trios Mucosal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) + (1|Litter) q<0.25")

make_taxa_dotplot(trios_het,
                  "Trios/final_taxonomy.tsv",
                  "Trios/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "Trios Mucosal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) + (1|Litter) q<0.25")


