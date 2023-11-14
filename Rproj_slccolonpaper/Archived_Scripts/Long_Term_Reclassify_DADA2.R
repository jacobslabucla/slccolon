
library(dada2)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/") # CHANGE to the directory containing the fastq files


# Assign taxonomy to ASVs using Silva database
seqtab <- read.delim("Long_Term/slc_LT_ASV_count_table_Oct2020.tsv", row.names = 1)
seqtab <- seqtab%>% select(-c("taxonomy"))
samples <- names(seqtab)
seqtab.nochim <- (t(seqtab))
row.names(seqtab.nochim) <-samples
seqtab.nochim <- as.data.frame(seqtab.nochim)
seqs <- names(seqtab.nochim)
class(seqs)

# Assign taxonomy to ASVs using Silva 138.2 database 

taxa <- assignTaxonomy(seqs, "C:/Users/Jacobs Laboratory/Desktop/16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/Jacobs Laboratory/Desktop/16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Export data so that it can be converted into BIOM
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
output<-cbind(t(seqtab.nochim), taxonomy)
#uniquesToFasta(seqtab.nochim, fout='C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/Long_Term/rep-seqs.fna', ids=colnames(seqtab.nochim))
write.table(output, "C:/Users/Jacobs Laboratory/Documents/JCYang/slccolonpaper/slccolon/Long_Term/SLT_ASV_table_Silva_v138_1.tsv", sep="\t", col.names=NA)

