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


setwd("/Users/kpark/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing/Porites_astreoides/DESeq/")  # The drag and drop from finder works in R, too.


#useful functions
#head() - prints out the top 6 lines
#dim() - prints the dimensions of a variable
#nrow() - returns the number of rows in a vector or matrix
# ?[functionName] - opens documentation describing the function

#read in your data to make counts table
countsTable <- read.delim('Mote_Past_Counts_full_cleannames.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

colSums(countsTable)
colSums(countsTable)[order(colSums(countsTable), decreasing=T)]

#make a table with conditions of each individual (e.g. "VA" and "RI")  There should be the same number of conditions described as there are samples in your data file, and in the same order.
#NOTE: It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.
#Here's an example
# Sample		Origin	SymState	Temp  Origin_SymState Origin_SymState_Temp
# RI_B_06_18		RI	B	18  RI_B  RI_B_18
# VA_B_07_14		VA	B	14  VA_B  VA_B_14
# VA_W_07_22		VA	W	22  VA_W  VA_W_22
# safest way to do this is in R

conds<-as.data.frame(matrix(unlist(strsplit(names(countsTable), "_")), ncol=4, byrow=T))
names(conds)<-c("Origin", "Geno","Temp","Timepoint")
conds$Sample<-names(countsTable)
conds$Geno<-factor(conds$Geno)
conds$Temp<-factor(conds$Temp)
conds$Timepoint<-factor(conds$Timepoint)
conds$Origin_Temp<-factor(paste(conds$Origin,conds$Temp, sep="_"))
conds$Origin_Temp_Timepoint<-factor(paste(conds$Origin,conds$Temp,conds$Timepoint, sep="_"))

# Samples Breakdown # 
table(conds$Origin)
table(conds$Geno)
table(conds$Temp)
table(conds$Timepoint)

##Check null counts per sample
prop.null <- apply(countsTable, 2, function(x) 100*mean(x==0))

barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="plots/OverallNullCountsv2.pdf",14, 14)
barplot(prop.null, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()

# how many genes have mean count >=3 
means = apply(countsTable, 1, mean)
table(means>=3)

# row sums for top 20 contigs
CountSums <- data.frame("ContigName"=row.names(countsTable), "CountSums"=rowSums(countsTable))
head(CountSums[order(CountSums$CountSums, decreasing = T),], n=20)

Top20Annos<-Annotations[Annotations$ContigName%in%row.names(head(CountSums[order(CountSums$CountSums, decreasing = T),], n=20)),]

# Calculating how many contigs have >= average 3 reads across all samples 
countsTable<-countsTable[means>=3,]
countsTable<-round(countsTable, digits=0)
dim(countsTable)


#Calculate proportion of contigs with zero reads within a sample 
prop.nullv2 <- apply(countsTable, 2, function(x) 100*mean(x==0))

prop.nullv3 <- apply(countsTable75, 2, function(x) 100*mean(x==0))
barplot(prop.nullv3, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1,
        xlim= c(0,80),
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')


#Calculate # of samples with > xx% null counts 
sum(prop.nullv2>90)
sum(prop.nullv2>80)
sum(prop.nullv2>75)
sum(prop.nullv2>70)



#create countsTable for 75%
countsTable75 <- countsTable[,prop.nullv2<75]

conds75 <- conds[conds$Sample%in%names(countsTable75),]

table(conds75$Origin)
table(conds75$Geno)
table(conds75$Temp)
table(conds75$Timepoint)

#make count data sets for 75%
dds75 <- DESeqDataSetFromMatrix(countData=countsTable75, colData=conds75, design=~ Origin_Temp_Timepoint)

dds75$Origin_Temp_Timepoint <- relevel(dds75$Origin_Temp_Timepoint, ref = "MoteIn_T0_T0")

dds75 <- DESeq(dds75)

dim(dds75) # dimensions of object

#transform counts to variance stabilized counts
vstCount75<-varianceStabilizingTransformation(dds75)

plotPCA(vstCount75, intgroup="Origin")


#transform counts to variance stablized counts
vstCounts<-varianceStabilizingTransformation(dds)


barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="plots/Mean3NullCountsv2.pdf",14, 14)
barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()


#make count data sets
dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=conds, design=~ Origin_Temp_Timepoint)

dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteIn_T0_T0")

dds <- DESeq(dds)

dim(dds) # dimmensions of object

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
colData(dds)

res <- results(dds)
head(res)
summary(res)

resInvOffT0 <- results(dds, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_0", "MoteInn_T0_0"))
head(resSymState)
summary(resSymState)


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
plotPCA(vstCounts, intgroup="Origin")
plotPCA(vstCounts, intgroup="Geno")
plotPCA(vstCounts, intgroup="Temp")

plotPCA(vstCounts, intgroup="Temp")

NamesToSubset<-conds[conds$Timepoint =="TF", "Sample"]
plotPCA(vstCounts[,colnames(vstCounts)%in%NamesToSubset], intgroup="Origin")

####PCAs for specific results sets####
i="TF"
j="T0"
Group<-conds[conds$Timepoint==i | conds$Timepoint==j,"Sample"]
pdf(file=paste0(i,"svs",j,"s_PCA.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()

i="TF"
Group<-conds[conds$Timepoint==i,"Sample"]
pdf(file=paste0("plots/", i,"s_PCA.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()

i="T5"
Group<-conds[conds$Timepoint==i,"Sample"]
#pdf(file=paste0("plots/",i,"s_PCA.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Temp")
dev.off()

i="T0"
Group<-conds[conds$Timepoint==i,"Sample"]
pdf(file=paste0("plots", i,"s_PCA.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()


i="T0"
j="MoteOff"
Group<-conds[conds$Temp==i | conds$Origin==j,"Sample"]
pdf(file=paste0(i,"svs",j,"s_PCA.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin")
dev.off()




#######################
######Jakes's Code#####
# vsd_Chlamy_N <- vst(dds_Chlamy_N, blind=FALSE)
#calculating the angle between two points in degrees

# Libraries 
library(ggplot2)
library(tidyverse)
library(cowplot)

# Create PCA df #
pcaData_Chlamy_N <- plotPCA(vstCounts, intgroup=c("Temp","Timepoint", "Origin"), returnData=TRUE)
percentVar_Chlamy_N <- round(100 * attr(pcaData_Chlamy_N, "percentVar"))
pcaData_Chlamy_N$Timepoint<-factor(pcaData_Chlamy_N$Timepoint, levels = c("TF","T0","T1","T2","T3","T4","T5"))
pcaData_Chlamy_N<-cbind(pcaData_Chlamy_N, pcaData_Chlamy_N[c(3:4,6:24,27,25:26,28:nrow(pcaData_Chlamy_N),1:2, 5), c(1,2,7)])
pcaData_Chlamy_N[c((nrow(pcaData_Chlamy_N)-2):nrow(pcaData_Chlamy_N),22),8:10]<-NA
colnames(pcaData_Chlamy_N)[8:10]<-c("PC1a", "PC2a", "name_a")

#calculating the angle between two points in degrees
pcaData_Chlamy_N$angle<-atan2( pcaData_Chlamy_N$PC2a-pcaData_Chlamy_N$PC2, pcaData_Chlamy_N$PC1a-pcaData_Chlamy_N$PC1) #*180/pi

# All Data With Arrows #

# Inshore
pcaData_MoteIn <- pcaData_Chlamy_N %>% 
  filter(Origin == "MoteIn")

a <- ggplot(pcaData_MoteIn, aes(PC1, PC2, fill=factor(Temp))) +
  geom_spoke(aes(angle = pcaData_MoteIn$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.95, shape = 21) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", "#FF61C3" )) +
  coord_fixed()+ggtitle("Inshore")+
  theme_light()+            
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.box.background = element_rect(), 
        legend.position = "none", legend.text = element_text(size = 8))

# Offshore
pcaData_MoteOff <- pcaData_Chlamy_N %>% 
  filter(Origin == "MoteOff")

b <- ggplot(pcaData_MoteOff, aes(PC1, PC2, fill=factor(Temp))) +
  geom_spoke(aes(angle = pcaData_MoteOff$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.95, shape = 22) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("Offshore")+
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", "#FF61C3" )) +
  theme_light()+            
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.box.background = element_rect(), 
        legend.position = "none", legend.text = element_text(size = 8))

# Combine Inshore and Offshore
plot_grid(a, b,  align = c("hv")) 

# Inshore and Offshore 
ggplot(pcaData_Chlamy_N, aes(PC1, PC2, fill=factor(Temp), shape=factor(Origin))) +
  geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=5,  alpha=.95) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  scale_shape_manual(values=c(21,22,23), name ="Nutrient")+
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("") +
  theme_light() +
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.box.background = element_rect(), 
        legend.position = "right", legend.text = element_text(size = 8))




# Remade PCA Plots #
pcaData<- pcaData_Chlamy_N %>%
  filter(Timepoint != "TF") # removing TF plot from figures 

# Inshore vs Offshore by Timepoint
ggplot(pcaData, aes(PC1, PC2, fill = factor(Temp), shape=factor(Origin))) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.85) + #ylim(-30,30)+xlim(-30,30)+
  #geom_text(aes(label = Temp), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF")) +
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+            
  facet_wrap(~Timepoint) +
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.position = "top")
ggsave("PCA_N.png", width = 8.5, height = 10, units = "in", dpi = "retina")

# Inshore vs Offshore by Temperature 
ggplot(pcaData, aes(PC1, PC2, fill = factor(Timepoint), shape=factor(Origin))) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=5,  alpha=.85) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  #scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF")) +
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+            
  facet_grid(Origin~Temp) +
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none")

# Attempted Arrrow Fig #
ggplot(pcaData_Chlamy_N, aes(PC1, PC2, color = factor(Timepoint))) +
  geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7) +
  geom_point(size=2,  alpha=.85) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2) +
  scale_shape_manual(values=c(21,22), name ="Origin")+
  #scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF")) +
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+   
  facet_wrap(~Temp) +
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14))



##############
# Individual Plots 
pcaData_T1<- pcaData_Chlamy_N %>%
  filter(Timepoint == "T1")

ggplot(pcaData_T1, aes(PC1, PC2, color = factor(Temp), shape=factor(Origin))) +
  geom_spoke(aes(angle = pcaData_T1$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7) +
  geom_point(size=7,  alpha=.85)

pcaData_30_In <- pcaData_Chlamy_N %>%
  filter(Temp == "30") %>%
  filter(Origin == "MoteIn")

ggplot(pcaData_30_In, aes(PC1, PC2, color = factor(Timepoint))) +
  geom_spoke(aes(angle = pcaData_30_In$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7) +
  geom_point(size=7,  alpha=.85)






######## 75% Plots ##########

# Create PCA df #
pcaData_75 <- plotPCA(vstCount75, intgroup=c("Temp","Timepoint", "Origin"), returnData=TRUE)
percentVar_75 <- round(100 * attr(pcaData_75, "percentVar"))
pcaData_75$Timepoint<-factor(pcaData_75$Timepoint, levels = c("TF","T0","T1","T2","T3","T4","T5"))
pcaData_75<-cbind(pcaData_75, pcaData_75[c(3:4,6:24,27,25:26,28:nrow(pcaData_75),1:2, 5), c(1,2,7)])
pcaData_75[c((nrow(pcaData_75)-2):nrow(pcaData_75),22),8:10]<-NA
colnames(pcaData_75)[8:10]<-c("PC1a", "PC2a", "name_a")

#calculating the angle between two points in degrees
pcaData_75$angle<-atan2( pcaData_75$PC2a-pcaData_75$PC2, pcaData_75$PC1a-pcaData_75$PC1) #*180/pi

pcaData_75 <- pcaData_75 %>%
  filter(Timepoint == "T0" | Timepoint == "TF")

ggplot(pcaData_75, aes(PC1, PC2, fill = factor(Temp), shape=factor(Origin))) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.85) +
  ylim(-50,60)+xlim(-50,150)+
  # geom_text(aes(label = Temp), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", 'pink')) +
  xlab(paste0("PC1: ",percentVar_75[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_75[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+      
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.position = "bottom")
ggsave("PCA_N.png", width = 8.5, height = 10, units = "in", dpi = "retina")


ggplot(pcaData_75, aes(PC1, PC2, fill = factor(Timepoint), shape=factor(Origin))) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=6,  alpha=.85) + #ylim(-30,30)+xlim(-30,30)+
  # geom_text(aes(label = Timepoint), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  #scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF")) +
  xlab(paste0("PC1: ",percentVar_75[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_75[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+            
  facet_grid(Origin~Temp) +
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none")

# Inshore
pcaData_75_MoteIn <- pcaData_75 %>% 
  filter(Origin == "MoteIn")

c <- ggplot(pcaData_75_MoteIn, aes(PC1, PC2, fill=factor(Temp))) +
  geom_spoke(aes(angle = pcaData_75_MoteIn$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.95, shape = 21) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  xlab(paste0("PC1: ",percentVar_75[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_75[2],"% variance")) + 
  xlim(c(-40,150)) +
  ylim(c(-60, 60)) +
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", "#FF61C3" )) +
  coord_fixed()+ggtitle("Inshore")+
  theme_light()+            
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.box.background = element_rect(), 
        legend.position = "none", legend.text = element_text(size = 8))

# Offshore
pcaData_75_MoteOff <- pcaData_75 %>% 
  filter(Origin == "MoteOff")

d <- ggplot(pcaData_75_MoteOff, aes(PC1, PC2, fill=factor(Temp))) +
  geom_spoke(aes(angle = pcaData_75_MoteOff$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
             radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.95, shape = 22) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  xlab(paste0("PC1: ",percentVar_75[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_75[2],"% variance")) + 
  xlim(c(-40,150)) +
  ylim(c(-60, 60)) +
  coord_fixed()+ggtitle("Offshore")+
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", "#FF61C3" )) +
  theme_light()+            
  theme(axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        title = element_text(size = 14, face = "italic"),
        legend.box.background = element_rect(), 
        legend.position = "none", legend.text = element_text(size = 8))


# Combine Inshore and Offshore
plot_grid(c, d,  align = c("hv")) 


pcaData_75 %>%
  filter(Timepoint == "TF")



















##############
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
Annotations<-read.delim("P_astreoides_Kenkel2013_totalannotated.txt")
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
