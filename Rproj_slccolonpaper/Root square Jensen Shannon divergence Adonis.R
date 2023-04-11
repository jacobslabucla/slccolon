library(vegan)
#setwd(choose.dir())
setwd("/home/julianne/Documents/SLC_project/Lumen_colononly/")
#data<-read.csv(choose.files(),header=T,row.names=1) # OTU biom filtered to 10-25% and unrarefied and converted to txt file then edited where column is sample and rows are OTU names similar to sPLS. Then save as cvs
#metadata<-read.csv(choose.files(),header=T,row.names=1)      #note: metadata file must be sorted to match the order of samples in the data file
#convert the L6.biom text into the CSV file 
data<- read.csv("ASV_Lumen_colononly.csv",header=T,row.names=1)
#data<-write.csv(x, file="Catia_ASV_count_table_no_singletons_L6.txt",append = FALSE,row.names=TRUE)
metadata<- read.csv("metadata_lumen_colononly.csv", header=T, row.names=1)

relativedata<-sweep(data,2,colSums(data),"/")

 dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
 }
## From: http://enterotype.embl.de/enterotypes.html

data.dist=dist.JSD(relativedata)

# to export distance matrix
x=data.dist
y=as.matrix(x)


write.table(y, "SLC_mucosa_SI_distance_matrix.txt", sep = "\t",row.names = TRUE, col.names = TRUE)
#distancematrixAGAIN <- read.table(file = 'Two Root_square_Jensen_Shannon_divergence_distance_matrix.txt', sep = " ", header = TRUE)
#write.csv(distancematrixAGAIN, file = "4thattempt_Two.csv")

# try exporting distance matrix as txt 
#write.table(y,"2ndattempt_Jan29_Root_square_Jensen_Shannon_divergence_distance_matrix.txt")

### To perform Adonis

data.adonis=adonis(data.dist ~ Litter + region + Genotype, data=metadata, permutations=10000)
data.adonis

##convert csv to txt file 
write.table(data, file = "Jan29_ReFormatted_Catia_ASV_count_table_no_singletons_L6.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)
write.table(data, file = "Jan29_ReFormatted_Sternini_Combined_Mapping.txt", sep = "\t",
            row.names = TRUE, col.names = TRUE)

#read in the distance matrix tsv file and convert to txt
distance_matrix <- read.table(file="Root_square_Jensen_Shannon_divergence_distance_matrix.txt", sep = '\t', header = TRUE)
#write.table(distance_matrix, file = "Jan29_Root_square_Jensen_Shannon_divergence_distance_matrix.txt", sep = "\t",
            #row.names = TRUE, col.names = TRUE)
distance_matrix_tabs <- write.csv(distance_matrix, file = "Juliannetab_Root_square_Jensen_Shannon_divergence_distance_matrix.txt", sep ="\t")



