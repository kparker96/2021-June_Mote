### This is an R script that implements the program DESeq2 for gene expression analysis.
### Much more information on the program and specific function (particularly for checking quality) 
### can be found here: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### The elements of this script were written by Melissa Pespeni and Daniel Barshis.

library(DESeq2)
library(gplots)
library(goseq)
library(GO.db)

setwd("/Users/dbarshis/dansstuff/Salmon_counts")  # The drag and drop from finder works in R, too.


####read in your data to make counts table####
countsTable <- read.delim('Mote_Past_Counts_full_cleannames.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

colSums(countsTable)
colSums(countsTable)[order(colSums(countsTable), decreasing=T)]

####make a table with conditions of each individual####
#(e.g. "VA" and "RI")  There should be the same number of conditions described as there are samples in your data file, and in the same order.
#NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.
#Here's an example
# Sample		Origin	SymState	Temp  Origin_SymState Origin_SymState_Temp
# RI_B_06_18		RI	B	18  RI_B  RI_B_18
# VA_B_07_14		VA	B	14  VA_B  VA_B_14
# VA_W_07_22		VA	W	22  VA_W  VA_W_22
# safest way to do this is in R

conds<-as.data.frame(matrix(unlist(strsplit(names(countsTable), "_")), ncol=4, byrow=T))
names(conds)<-c("Origin","Geno","Temp","Timepoint")
conds$Sample<-names(countsTable)
conds$Geno<-factor(conds$Geno)
conds$Temp<-factor(conds$Temp)
conds$Timepoint<-factor(conds$Timepoint)
conds$Origin_Temp<-factor(paste(conds$Origin,conds$Temp, sep="_"))
conds$Origin_Temp_Timepoint<-factor(paste(conds$Origin,conds$Temp,conds$Timepoint, sep="_"))

####Check replications####
table(conds$Origin)
table(conds$Geno)
table(conds$Temp)
table(conds$Timepoint)

####Check null counts per sample####
prop.null <- apply(countsTable, 2, function(x) 100*mean(x==0))

barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="OverallNullCountsv2.pdf",14, 14)
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()

#### how many genes have mean count >=3 ####
means = apply(countsTable, 1, mean)
table(means>=3)

#### row sums for top 20 contigs####
CountSums <- data.frame("ContigName"=row.names(countsTable), "CountSums"=rowSums(countsTable))
head(CountSums[order(CountSums$CountSums, decreasing = T),], n=20)

#Top20Annos<-Annotations[Annotations$ContigName%in%row.names(head(CountSums[order(CountSums$CountSums, decreasing = T),], n=20)),]


countsTable<-countsTable[means>=3,]
countsTable<-round(countsTable, digits=0)
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

####make count data sets####
dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=conds, design=~ Origin_Temp_Timepoint)

#### Specific contrasts ####
#### MoteOff_T0_T0 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_T0_T0")

ddsOffT0 <- DESeq(dds)

dim(ddsOffT0)
resultsNames(ddsOffT0)

resInvOffT0 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteIn_T0_T0"))
write.table(resInvOffT0, file="resInvOffT0.txt", quote=F, row.names=T, sep="\t")
head(resInvOffT0)
summary(resInvOffT0)

resOffT0vOff_TF_TF <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_TF_TF"))
write.table(resOffT0vOff_TF_TF, file="resOffT0vOff_TF_TF.txt", quote=F, row.names=T, sep="\t")
resOffT0vOff_30_T1 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_30_T1"))
write.table(resOffT0vOff_30_T1, file="resOffT0vOff_30_T1.txt", quote=F, row.names=T, sep="\t")
resOffT0vOff_30_T2 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_30_T2"))
write.table(resOffT0vOff_30_T2, file="resOffT0vOff_30_T2.txt", quote=F, row.names=T, sep="\t")
resOffT0vOff_30_T3 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_30_T3"))
write.table(resOffT0vOff_30_T3, file="resOffT0vOff_30_T3.txt", quote=F, row.names=T, sep="\t")
resOffT0vOff_30_T4 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_30_T4"))
write.table(resOffT0vOff_30_T4, file="resOffT0vOff_30_T4.txt", quote=F, row.names=T, sep="\t")
resOffT0vOff_30_T5 <- results(ddsOffT0, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteOff_30_T5"))
write.table(resOffT0vOff_30_T5, file="resOffT0vOff_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resInvOffT0)
summary(resOffT0vOff_TF_TF)
summary(resOffT0vOff_30_T1)
summary(resOffT0vOff_30_T2)
summary(resOffT0vOff_30_T3)
summary(resOffT0vOff_30_T4)
summary(resOffT0vOff_30_T5)
date()
#### MoteIn_T0_T0 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_T0_T0")

ddsInT0 <- DESeq(dds)

dim(ddsInT0)
resultsNames(ddsInT0)

resInT0vIn_TF_TF <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_TF_TF"))
write.table(resInT0vIn_TF_TF, file="resInT0vIn_TF_TF.txt", quote=F, row.names=T, sep="\t")
resInT0vIn_30_T1 <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_30_T1"))
write.table(resInT0vIn_30_T1, file="resInT0vIn_30_T1.txt", quote=F, row.names=T, sep="\t")
resInT0vIn_30_T2 <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_30_T2"))
write.table(resInT0vIn_30_T2, file="resInT0vIn_30_T2.txt", quote=F, row.names=T, sep="\t")
resInT0vIn_30_T3 <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_30_T3"))
write.table(resInT0vIn_30_T3, file="resInT0vIn_30_T3.txt", quote=F, row.names=T, sep="\t")
resInT0vIn_30_T4 <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_30_T4"))
write.table(resInT0vIn_30_T4, file="resInT0vIn_30_T4.txt", quote=F, row.names=T, sep="\t")
resInT0vIn_30_T5 <- results(ddsInT0, contrast=c("Origin_Temp_Timepoint", "MoteIn_T0_T0", "MoteIn_30_T5"))
write.table(resInT0vIn_30_T5, file="resInT0vIn_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resInT0vIn_TF_TF)
summary(resInT0vIn_30_T1)
summary(resInT0vIn_30_T2)
summary(resInT0vIn_30_T3)
summary(resInT0vIn_30_T4)
summary(resInT0vIn_30_T5)
date()

#### MoteOff_30_T1 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_30_T1")

ddsOff_30_T1 <- DESeq(dds)

dim(ddsOff_30_T1)
resultsNames(ddsOff_30_T1)

resOff_30_T1vIn_30_T1 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteIn_30_T1"))
write.table(resOff_30_T1vIn_30_T1, file="resOff_30_T1vIn_30_T1.txt", quote=F, row.names=T, sep="\t")
head(resOff_30_T1vIn_30_T1)
summary(resOff_30_T1vIn_30_T1)

resOff_30_T1vOff_34_T1 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_34_T1"))
write.table(resOff_30_T1vOff_34_T1, file="resOff_30_T1vOff_34_T1.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_37_T1 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_37_T1"))
write.table(resOff_30_T1vOff_37_T1, file="resOff_30_T1vOff_37_T1.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_39_T1 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_39_T1"))
write.table(resOff_30_T1vOff_39_T1, file="resOff_30_T1vOff_39_T1.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_30_T2 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_30_T2"))
write.table(resOff_30_T1vOff_30_T2, file="resOff_30_T1vOff_30_T2.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_30_T3 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_30_T3"))
write.table(resOff_30_T1vOff_30_T3, file="resOff_30_T1vOff_30_T3.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_30_T4 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_30_T4"))
write.table(resOff_30_T1vOff_30_T4, file="resOff_30_T1vOff_30_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T1vOff_30_T5 <- results(ddsOff_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T1", "MoteOff_30_T5"))
write.table(resOff_30_T1vOff_30_T5, file="resOff_30_T1vOff_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resOff_30_T1vIn_30_T1)
summary(resOff_30_T1vOff_34_T1)
summary(resOff_30_T1vOff_37_T1)
summary(resOff_30_T1vOff_39_T1)
summary(resOff_30_T1vOff_30_T2)
summary(resOff_30_T1vOff_30_T3)
summary(resOff_30_T1vOff_30_T4)
summary(resOff_30_T1vOff_30_T5)
date()

#### MoteIn_30_T1 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_30_T1")

ddsIn_30_T1 <- DESeq(dds)

dim(ddsIn_30_T1)
resultsNames(ddsIn_30_T1)

resIn_30_T1vIn_34_T1 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_34_T1"))
write.table(resIn_30_T1vIn_34_T1, file="resIn_30_T1vIn_34_T1.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_37_T1 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_37_T1"))
write.table(resIn_30_T1vIn_37_T1, file="resIn_30_T1vIn_37_T1.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_39_T1 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_39_T1"))
write.table(resIn_30_T1vIn_39_T1, file="resIn_30_T1vIn_39_T1.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_30_T2 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_30_T2"))
write.table(resIn_30_T1vIn_30_T2, file="resIn_30_T1vIn_30_T2.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_30_T3 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_30_T3"))
write.table(resIn_30_T1vIn_30_T3, file="resIn_30_T1vIn_30_T3.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_30_T4 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_30_T4"))
write.table(resIn_30_T1vIn_30_T4, file="resIn_30_T1vIn_30_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T1vIn_30_T5 <- results(ddsIn_30_T1, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T1", "MoteIn_30_T5"))
write.table(resIn_30_T1vIn_30_T5, file="resIn_30_T1vIn_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resIn_30_T1vIn_34_T1)
summary(resIn_30_T1vIn_37_T1)
summary(resIn_30_T1vIn_39_T1)
summary(resIn_30_T1vIn_30_T2)
summary(resIn_30_T1vIn_30_T3)
summary(resIn_30_T1vIn_30_T4)
summary(resIn_30_T1vIn_30_T5)
date()

#### MoteOff_30_T2 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_30_T2")

ddsOff_30_T2 <- DESeq(dds)

dim(ddsOff_30_T2)
resultsNames(ddsOff_30_T2)

resOff_30_T2vIn_30_T2 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteIn_30_T2"))
write.table(resOff_30_T2vIn_30_T2, file="resOff_30_T2vIn_30_T2.txt", quote=F, row.names=T, sep="\t")
head(resOff_30_T2vIn_30_T2)
summary(resOff_30_T2vIn_30_T2)

resOff_30_T2vOff_34_T2 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_34_T2"))
write.table(resOff_30_T2vOff_34_T2, file="resOff_30_T2vOff_34_T2.txt", quote=F, row.names=T, sep="\t")
resOff_30_T2vOff_37_T2 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_37_T2"))
write.table(resOff_30_T2vOff_37_T2, file="resOff_30_T2vOff_37_T2.txt", quote=F, row.names=T, sep="\t")
resOff_30_T2vOff_39_T2 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_39_T2"))
write.table(resOff_30_T2vOff_39_T2, file="resOff_30_T2vOff_39_T2.txt", quote=F, row.names=T, sep="\t")
resOff_30_T2vOff_30_T3 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_30_T3"))
write.table(resOff_30_T2vOff_30_T3, file="resOff_30_T2vOff_30_T3.txt", quote=F, row.names=T, sep="\t")
resOff_30_T2vOff_30_T4 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_30_T4"))
write.table(resOff_30_T2vOff_30_T4, file="resOff_30_T2vOff_30_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T2vOff_30_T5 <- results(ddsOff_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T2", "MoteOff_30_T5"))
write.table(resOff_30_T2vOff_30_T5, file="resOff_30_T2vOff_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resOff_30_T2vIn_30_T2)
summary(resOff_30_T2vOff_34_T2)
summary(resOff_30_T2vOff_37_T2)
summary(resOff_30_T2vOff_39_T2)
summary(resOff_30_T2vOff_30_T3)
summary(resOff_30_T2vOff_30_T4)
summary(resOff_30_T2vOff_30_T5)
date()

#### MoteIn_30_T2 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_30_T2")

ddsIn_30_T2 <- DESeq(dds)

dim(ddsIn_30_T2)
resultsNames(ddsIn_30_T2)

resIn_30_T2vIn_34_T2 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_34_T2"))
write.table(resIn_30_T2vIn_34_T2, file="resIn_30_T2vIn_34_T2.txt", quote=F, row.names=T, sep="\t")
resIn_30_T2vIn_37_T2 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_37_T2"))
write.table(resIn_30_T2vIn_37_T2, file="resIn_30_T2vIn_37_T2.txt", quote=F, row.names=T, sep="\t")
resIn_30_T2vIn_39_T2 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_39_T2"))
write.table(resIn_30_T2vIn_39_T2, file="resIn_30_T2vIn_39_T2.txt", quote=F, row.names=T, sep="\t")
resIn_30_T2vIn_30_T3 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_30_T3"))
write.table(resIn_30_T2vIn_30_T3, file="resIn_30_T2vIn_30_T3.txt", quote=F, row.names=T, sep="\t")
resIn_30_T2vIn_30_T4 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_30_T4"))
write.table(resIn_30_T2vIn_30_T4, file="resIn_30_T2vIn_30_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T2vIn_30_T5 <- results(ddsIn_30_T2, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T2", "MoteIn_30_T5"))
write.table(resIn_30_T2vIn_30_T5, file="resIn_30_T2vIn_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resIn_30_T2vIn_34_T2)
summary(resIn_30_T2vIn_37_T2)
summary(resIn_30_T2vIn_39_T2)
summary(resIn_30_T2vIn_30_T3)
summary(resIn_30_T2vIn_30_T4)
summary(resIn_30_T2vIn_30_T5)

date()

#### MoteOff_30_T3 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_30_T3")

ddsOff_30_T3 <- DESeq(dds)

dim(ddsOff_30_T3)
resultsNames(ddsOff_30_T3)

resOff_30_T3vIn_30_T3 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteIn_30_T3"))
write.table(resOff_30_T3vIn_30_T3, file="resOff_30_T3vIn_30_T3.txt", quote=F, row.names=T, sep="\t")
head(resOff_30_T3vIn_30_T3)
summary(resOff_30_T3vIn_30_T3)

resOff_30_T3vOff_34_T3 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteOff_34_T3"))
write.table(resOff_30_T3vOff_34_T3, file="resOff_30_T3vOff_34_T3.txt", quote=F, row.names=T, sep="\t")
resOff_30_T3vOff_37_T3 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteOff_37_T3"))
write.table(resOff_30_T3vOff_37_T3, file="resOff_30_T3vOff_37_T3.txt", quote=F, row.names=T, sep="\t")
resOff_30_T3vOff_39_T3 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteOff_39_T3"))
write.table(resOff_30_T3vOff_39_T3, file="resOff_30_T3vOff_39_T3.txt", quote=F, row.names=T, sep="\t")
resOff_30_T3vOff_30_T4 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteOff_30_T4"))
write.table(resOff_30_T3vOff_30_T4, file="resOff_30_T3vOff_30_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T3vOff_30_T5 <- results(ddsOff_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T3", "MoteOff_30_T5"))
write.table(resOff_30_T3vOff_30_T5, file="resOff_30_T3vOff_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resOff_30_T3vOff_34_T3)
summary(resOff_30_T3vOff_37_T3)
summary(resOff_30_T3vOff_39_T3)
summary(resOff_30_T3vOff_30_T4)
summary(resOff_30_T3vOff_30_T5)
date()

#### MoteIn_30_T3 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_30_T3")

ddsIn_30_T3 <- DESeq(dds)

dim(ddsIn_30_T3)
resultsNames(ddsIn_30_T3)

resIn_30_T3vIn_34_T3 <- results(ddsIn_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T3", "MoteIn_34_T3"))
write.table(resIn_30_T3vIn_34_T3, file="resIn_30_T3vIn_34_T3.txt", quote=F, row.names=T, sep="\t")
resIn_30_T3vIn_37_T3 <- results(ddsIn_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T3", "MoteIn_37_T3"))
write.table(resIn_30_T3vIn_37_T3, file="resIn_30_T3vIn_37_T3.txt", quote=F, row.names=T, sep="\t")
resIn_30_T3vIn_39_T3 <- results(ddsIn_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T3", "MoteIn_39_T3"))
write.table(resIn_30_T3vIn_39_T3, file="resIn_30_T3vIn_39_T3.txt", quote=F, row.names=T, sep="\t")
resIn_30_T3vIn_30_T4 <- results(ddsIn_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T3", "MoteIn_30_T4"))
write.table(resIn_30_T3vIn_30_T4, file="resIn_30_T3vIn_30_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T3vIn_30_T5 <- results(ddsIn_30_T3, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T3", "MoteIn_30_T5"))
write.table(resIn_30_T3vIn_30_T5, file="resIn_30_T3vIn_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resIn_30_T3vIn_34_T3)
summary(resIn_30_T3vIn_37_T3)
summary(resIn_30_T3vIn_39_T3)
summary(resIn_30_T3vIn_30_T4)
summary(resIn_30_T3vIn_30_T5)

date()

#### MoteOff_30_T4 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_30_T4")

ddsOff_30_T4 <- DESeq(dds)

dim(ddsOff_30_T4)
resultsNames(ddsOff_30_T4)

resOff_30_T4vIn_30_T4 <- results(ddsOff_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T4", "MoteIn_30_T4"))
write.table(resOff_30_T4vIn_30_T4, file="resOff_30_T4vIn_30_T4.txt", quote=F, row.names=T, sep="\t")
head(resOff_30_T4vIn_30_T4)
summary(resOff_30_T4vIn_30_T4)

resOff_30_T4vOff_34_T4 <- results(ddsOff_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T4", "MoteOff_34_T4"))
write.table(resOff_30_T4vOff_34_T4, file="resOff_30_T4vOff_34_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T4vOff_37_T4 <- results(ddsOff_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T4", "MoteOff_37_T4"))
write.table(resOff_30_T4vOff_37_T4, file="resOff_30_T4vOff_37_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T4vOff_39_T4 <- results(ddsOff_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T4", "MoteOff_39_T4"))
write.table(resOff_30_T4vOff_39_T4, file="resOff_30_T4vOff_39_T4.txt", quote=F, row.names=T, sep="\t")
resOff_30_T4vOff_30_T5 <- results(ddsOff_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T4", "MoteOff_30_T5"))
write.table(resOff_30_T4vOff_30_T5, file="resOff_30_T4vOff_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resOff_30_T4vOff_34_T4)
summary(resOff_30_T4vOff_37_T4)
summary(resOff_30_T4vOff_39_T4)
summary(resOff_30_T4vOff_30_T5)
date()

#### MoteIn_30_T4 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_30_T4")

ddsIn_30_T4 <- DESeq(dds)

dim(ddsIn_30_T4)
resultsNames(ddsIn_30_T4)

resIn_30_T4vIn_34_T4 <- results(ddsIn_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T4", "MoteIn_34_T4"))
write.table(resIn_30_T4vIn_34_T4, file="resIn_30_T4vIn_34_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T4vIn_37_T4 <- results(ddsIn_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T4", "MoteIn_37_T4"))
write.table(resIn_30_T4vIn_37_T4, file="resIn_30_T4vIn_37_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T4vIn_39_T4 <- results(ddsIn_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T4", "MoteIn_39_T4"))
write.table(resIn_30_T4vIn_39_T4, file="resIn_30_T4vIn_39_T4.txt", quote=F, row.names=T, sep="\t")
resIn_30_T4vIn_30_T5 <- results(ddsIn_30_T4, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T4", "MoteIn_30_T5"))
write.table(resIn_30_T4vIn_30_T5, file="resIn_30_T4vIn_30_T5.txt", quote=F, row.names=T, sep="\t")

summary(resIn_30_T4vIn_34_T4)
summary(resIn_30_T4vIn_37_T4)
summary(resIn_30_T4vIn_39_T4)
summary(resIn_30_T4vIn_30_T5)

date()

#### MoteOff_30_T5 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_30_T5")

ddsOff_30_T5 <- DESeq(dds)

dim(ddsOff_30_T5)
resultsNames(ddsOff_30_T5)

resOff_30_T5vIn_30_T5 <- results(ddsOff_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T5", "MoteIn_30_T5"))
write.table(resOff_30_T5vIn_30_T5, file="resOff_30_T5vIn_30_T5.txt", quote=F, row.names=T, sep="\t")
head(resOff_30_T5vIn_30_T5)
summary(resOff_30_T5vIn_30_T5)

resOff_30_T5vOff_34_T5 <- results(ddsOff_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T5", "MoteOff_34_T5"))
write.table(resOff_30_T5vOff_34_T5, file="resOff_30_T5vOff_34_T5.txt", quote=F, row.names=T, sep="\t")
resOff_30_T5vOff_37_T5 <- results(ddsOff_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T5", "MoteOff_37_T5"))
write.table(resOff_30_T5vOff_37_T5, file="resOff_30_T5vOff_37_T5.txt", quote=F, row.names=T, sep="\t")
resOff_30_T5vOff_39_T5 <- results(ddsOff_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteOff_30_T5", "MoteOff_39_T5"))
write.table(resOff_30_T5vOff_39_T5, file="resOff_30_T5vOff_39_T5.txt", quote=F, row.names=T, sep="\t")

summary(resOff_30_T5vOff_34_T5)
summary(resOff_30_T5vOff_37_T5)
summary(resOff_30_T5vOff_39_T5)
date()

#### MoteIn_30_T5 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_30_T5")

ddsIn_30_T5 <- DESeq(dds)

dim(ddsIn_30_T5)
resultsNames(ddsIn_30_T5)

resIn_30_T5vIn_34_T5 <- results(ddsIn_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T5", "MoteIn_34_T5"))
write.table(resIn_30_T5vIn_34_T5, file="resIn_30_T5vIn_34_T5.txt", quote=F, row.names=T, sep="\t")
resIn_30_T5vIn_37_T5 <- results(ddsIn_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T5", "MoteIn_37_T5"))
write.table(resIn_30_T5vIn_37_T5, file="resIn_30_T5vIn_37_T5.txt", quote=F, row.names=T, sep="\t")
resIn_30_T5vIn_39_T5 <- results(ddsIn_30_T5, contrast=c("Origin_Temp_Timepoint", "MoteIn_30_T5", "MoteIn_39_T5"))
write.table(resIn_30_T5vIn_39_T5, file="resIn_30_T5vIn_39_T5.txt", quote=F, row.names=T, sep="\t")

summary(resIn_30_T5vIn_34_T5)
summary(resIn_30_T5vIn_37_T5)
summary(resIn_30_T5vIn_39_T5)
date()

####transform counts to variance stablized counts####
vstCounts<-varianceStabilizingTransformation(dds)
