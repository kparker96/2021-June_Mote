### This is an R script that implements the program DESeq2 for gene expression analysis.
### Much more information on the program and specific function (particularly for checking quality) 
### can be found here: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### The elements of this script were written by Melissa Pespeni and Daniel Barshis.

#Only need to do this the first time to install the package
#for R version 4.0
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("goseq")

library(DESeq2)
library(gplots)
library(goseq)
library(GO.db)


setwd("YOURWORKINGDIRECTORY")  # The drag and drop from finder works in R, too.


#useful functions
#head() - prints out the top 6 lines
#dim() - prints the dimensions of a variable
#nrow() - returns the number of rows in a vector or matrix
# ?[functionName] - opens documentation describing the function

#read in your data to make counts table
countsTable <- read.delim('djbFullCounts_summed.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

colSums(countsTable)
colSums(countsTable)[order(colSums(countsTable))]

table(substring(names(countsTable),1,2))
table(substring(names(countsTable),4,4))
table(substring(names(countsTable),9,10))
table(paste(substring(names(countsTable),1,2),substring(names(countsTable),4,4), sep="_"))
table(paste(substring(names(countsTable),1,2),substring(names(countsTable),4,4), substring(names(countsTable),9,10), sep="_"))

LowCounts<-c("VA_W_02_22","VA_W_01_22","RI_B_03_18", "RI_W_01_14", "RI_W_04_22", "VA_W_03_18")
names(countsTable)%in%LowCounts
!(names(countsTable)%in%LowCounts)
countsTable<-countsTable[,!(names(countsTable)%in%LowCounts)]
dim(countsTable)

#make a table with conditions of each individual (e.g. "VA" and "RI")  There should be the same number of conditions described as there are samples in your data file, and in the same order.
#NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.
#Here's an example
# Sample		Origin	SymState	Temp  Origin_SymState Origin_SymState_Temp
# RI_B_06_18		RI	B	18  RI_B  RI_B_18
# VA_B_07_14		VA	B	14  VA_B  VA_B_14
# VA_W_07_22		VA	W	22  VA_W  VA_W_22
# safest way to do this is in R

conds<-data.frame("Sample"=names(countsTable))
conds$Origin<-factor(substring(conds$Sample,1,2))
conds$SymState<-factor(substring(conds$Sample,4,4))
conds$Temp<-factor(substring(conds$Sample,9,10))
conds$Origin_SymState<-factor(paste(conds$Origin,conds$SymState, sep="_"))
conds$Origin_SymState_Temp<-factor(paste(conds$Origin,conds$SymState,conds$Temp, sep="_"))

##Check null counts per sample
prop.null <- apply(countsTable, 2, function(x) 100*mean(x==0))

barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="OverallNullCountsv2.pdf",14, 14)
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()

# how many genes have mean count >=3 
means = apply(countsTable, 1, mean)
table(means>=3)

countsTable<-countsTable[means>=3,]
dim(countsTable)

prop.nullv2 <- apply(countsTable, 2, function(x) 100*mean(x==0))

barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="Mean3NullCountsv2.pdf",14, 14)
barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()

#make count data sets
dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=conds, design=~ Origin + SymState + Temp + Origin:SymState)

dds <- DESeq(dds)

dim(dds)

#transform counts to variance stablized counts
vstCounts<-varianceStabilizingTransformation(dds)

#plot cluster dendrogram across samples
dists<-dist(t(assay(vstCounts)))
plot(hclust(dists))
heatmap.2(as.matrix(dists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10))

###Figure out which contrast you want to examine (i.e. which two groups do you want to compare)
resultsNames(dds)

res <- results(dds)
head(res)
summary(res)

resSymState <- results(dds, contrast=c("SymState", "B", "W"))
head(resSymState)
summary(resSymState)

resOrigin <- results(dds, contrast=c("Origin", "VA", "RI"))
head(resOrigin)
summary(resOrigin)

res14v18 <- results(dds, contrast=c("Temp", "14", "18"))
head(res14v18)
summary(res14v18)

res14v22 <- results(dds, contrast=c("Temp", "14", "22"))
head(res14v22)
summary(res14v22)

res18v22 <- results(dds, contrast=c("Temp", "18", "22"))
head(res18v22)
summary(res18v22)

#count the number of significantly differentially expressed genes
sum(resOrigin$padj < 0.1, na.rm =T)
sum(resOrigin$pvalue < 0.05, na.rm =T)
											 
# Make a counts table that is scaled by the size factors
scaledcounts = counts(dds, normalize=T)
head(scaledcounts)

#building heat map data
head(scaledcounts)
genes4heatmap<-resOrigin[resOrigin$pvalue <0.05 & !is.na(resOrigin$pvalue),]
names(genes4heatmap)
head(genes4heatmap)
dim(genes4heatmap)
data4heatmap<-scaledcounts[row.names(scaledcounts)%in%row.names(genes4heatmap),]
dim(data4heatmap)
head(data4heatmap)

temp = as.matrix(rowMeans(data4heatmap))
head(temp)
scaledmatrix<-matrix(data=temp, nrow=nrow(data4heatmap), ncol=ncol(data4heatmap), byrow=FALSE)
data4heatmapscaled = data4heatmap/scaledmatrix
head(data4heatmapscaled)

dim(data4heatmapscaled)

pairs.breaks <- seq(0, 3.0, by=0.1)
length(pairs.breaks)
mycol <- colorpanel(n=30, low="black", high="yellow") 

pdf(file="XXXXXX_byrow.pdf",7,7)
heatmap.2(data.matrix(data4heatmapscaled), Rowv=T, Colv=F, dendrogram = c("row"), scale="none", keysize=1, breaks=pairs.breaks, col=mycol, trace = "none", symkey = F, density.info = "density", colsep=c(24), sepcolor=c("white"), sepwidth=c(.1,.1), margins=c(10,10), labRow=F)
dev.off()

pdf(file="XXXXXX_bycolumn.pdf",7,7)
heatmap.2(data.matrix(data4heatmapscaled), Rowv=T, Colv=T, dendrogram = c("col"), scale="none", keysize=1, breaks=pairs.breaks, col=mycol, trace = "none", symkey = F, density.info = "density", colsep=c(24), sepcolor=c("white"), sepwidth=c(.1,.1), margins=c(10,10), labRow=F)
dev.off()

###PCAs for specific results sets####
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="Origin")
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="SymState")
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="Temp")
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="Origin_SymState")
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="Origin_SymState_Temp")

pdf(file="866_RIvsVA_PCAv2.pdf",14, 14)
plotPCA(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),], intgroup="Origin")
dev.off()

prop.nullv3 <- apply(countsTable[rownames(countsTable)%in%row.names(genes4heatmap),], 2, function(x) 100*mean(x==0))

barplot(prop.nullv3, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="ResNullCountsv2.pdf",14, 14)
barplot(prop.nullv3, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()

distsv2<-dist(t(assay(vstCounts[rownames(vstCounts)%in%row.names(genes4heatmap),])))
plot(hclust(distsv2))
heatmap.2(as.matrix(distsv2), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10))


#### GOSeq analysis ####
res<-resOrigin
Annotations<-read.delim("Assembly_GoodCoralSymbiont_suffixed_totalannotated.txt")
head(res)
sum(res$log2FoldChange!=0)

genes<-as.integer(p.adjust(res$pvalue[res$log2FoldChange!=0 & !is.na(res$pvalue)], method="BH")<.05)
names(genes)<-row.names(res[res$log2FoldChange!=0 & !is.na(res$pvalue),])
table(genes)

genes<-as.integer(res$pvalue[res$log2FoldChange!=0 & !is.na(res$pvalue)]<.05)
names(genes)<-row.names(res[res$log2FoldChange!=0 & !is.na(res$pvalue),])
table(genes)

goterms<-strsplit(as.character(Annotations$GO), split=" // ")
names(goterms)<-Annotations$ContigName

pwf<-nullp(genes, bias.data=Annotations[Annotations$ContigName%in%names(genes),"ContigLength"])
GOGOGO<-goseq(pwf,gene2cat=goterms)
GOGOGO$padj<-p.adjust(GOGOGO$over_represented_pvalue, method="fdr")

sum(GOGOGO$padj<0.2, na.rm=T)
sum(GOGOGO$padj<0.1, na.rm=T)
sum(GOGOGO$padj<0.05, na.rm=T)
GOGOGO[GOGOGO$padj<0.05,]

GOTERM[["SOMEGOTERM"]]

EnGOS<-GOGOGO$category[p.adjust(GOGOGO$over_represented_pvalue, method="BH")<.05]

for(go in EnGOS){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")}
