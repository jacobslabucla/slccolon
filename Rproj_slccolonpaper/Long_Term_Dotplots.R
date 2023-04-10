
library(dada2)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files


## Make a taxa dotplot for SLT q<0.25 --
data<-read.table("Long_Term/Luminal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")

# Luminal Colon 
make_taxa_dotplot(data,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.25", value)

make_taxa_dotplot(trios_mut,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.25", Phylum)

make_taxa_dotplot(trios_het,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.25",Phylum)

# Mucosal Colon
data<-read.table("Long_Term/Mucosal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.25)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")


make_taxa_dotplot(data,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon: ASV ~ Site + Sex + Genotype +(1|MouseID) q<0.25", value)

make_taxa_dotplot(trios_mut,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.25", Phylum)

make_taxa_dotplot(trios_het,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.25", Phylum)

## Make a taxa dotplot for SLT q<0.05 --
data<-read.table("Long_Term/Luminal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")

# Luminal Colon 
make_taxa_dotplot(data,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.05")

make_taxa_dotplot(trios_mut,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.05")

make_taxa_dotplot(trios_het,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Luminal Colon-ASV.RDS",
                  "SLT Luminal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.05")

# Mucosal Colon
data<-read.table("Long_Term/Mucosal Colon_ASV_Maaslin2_Site_Sex_Genotype_1-MouseID/significant_results.tsv", header=TRUE)
data <- data %>% filter(qval <0.05)
data <- data %>% filter(metadata=="Genotype")

trios_mut <- data %>% filter(value=="MUT")
trios_het <- data %>% filter(value=="HET")

make_taxa_dotplot(data,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon: ASV ~ Site + Sex + Genotype +(1|MouseID) q<0.05", value)

make_taxa_dotplot(trios_mut,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon MUT: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.05", Phylum)

make_taxa_dotplot(trios_het,
                  "Long_Term/final_taxonomy.tsv",
                  "Long_Term/Relative_Abundance-Mucosal Colon-ASV.RDS",
                  "SLT Mucosal Colon HET: ASV ~ Site + Sex + Genotype + (1|MouseID) q<0.05", Phylum)



