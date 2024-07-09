### This is an R script that implements the program DESeq2 for gene expression analysis.
### Much more information on the program and specific function (particularly for checking quality) 
### can be found here: https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### The elements of this script were written by Melissa Pespeni and Daniel Barshis.

#Only need to do this the first time to install the package
#for R version 4.0
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
#BiocManager::install("DESeq2")
#BiocManager::install("geneLenDataBase")
#BiocManager::install("goseq")

library(DESeq2) # install.packages("DESeq2")
library(gplots)
library(goseq) # install.packages("goseq")
library(GO.db) # install.packages("GO.db")
library(ggplot2)
library(dplyr)

# setwd("/Users/kpark/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing/Porites_astreoides/Pilot_Fail_New_Concat_Compiled/")  # The drag and drop from finder works in R, too.

#useful functions
#head() - prints out the top 6 lines
#dim() - prints the dimensions of a variable
#nrow() - returns the number of rows in a vector or matrix
# ?[functionName] - opens documentation describing the function

#read in your data to make counts table
countsTable <- read.delim('data/Sequencing/Porites_astreoides/Pilot_Fail_New_Concat_Compiled/Dedup_Mote_Past_Counts_full_cleannames.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1) 

countsTable <- dplyr::select(countsTable, -read_sum, -read_mean)

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

Annotations<-read.delim("data/Sequencing/Porites_astreoides/Pilot_Fail_New_Concat_Compiled/P_astreoides_Kenkel2013_totalannotated.txt")

Top20Annos<-Annotations[Annotations$ContigName%in%row.names(head(CountSums[order(CountSums$CountSums, decreasing = T),], n=20)),]

# Calculating how many contigs have >= average 3 reads across all samples 
countsTable<-countsTable[means>=3,]
countsTable<-round(countsTable, digits=0)
dim(countsTable)


#Calculate proportion of contigs with zero reads within a sample 
prop.nullv2 <- apply(countsTable, 2, function(x) 100*mean(x==0))


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

prop.nullv3 <- apply(countsTable75, 2, function(x) 100*mean(x==0))
barplot(prop.nullv3, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1,
        xlim= c(0,80),
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')


#make count data sets for 75%
dds75 <- DESeq2::DESeqDataSetFromMatrix(countData=countsTable75, colData=conds75, design=~ Origin_Temp_Timepoint)

dds75$Origin_Temp_Timepoint <- relevel(dds75$Origin_Temp_Timepoint, ref = "MoteIn_T0_T0")

dds75 <- DESeq(dds75)

dim(dds75) # dimensions of object

#transform counts to variance stabilized counts
vstCount75<-varianceStabilizingTransformation(dds75)

plotPCA(vstCount75, intgroup="Origin")

barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')

pdf(file="plots/Mean3NullCountsv2.pdf",14, 14)
barplot(prop.nullv2, main="Percentage of null counts per sample", 
        horiz=TRUE, cex.names=0.5, las=1, 
        col=conds$Origin_SymState_Temp, ylab='Samples', xlab='% of null counts')
dev.off()


#make count data sets
dds <- DESeq2::DESeqDataSetFromMatrix(countData=countsTable, colData=conds, design=~ Origin_Temp_Timepoint)

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

resInvOffT0 <- results(dds, contrast=c("Origin_Temp_Timepoint", "MoteOff_T0_T0", "MoteIn_T0_T0"))
head(resInvOffT0)
summary(resInvOffT0)

#count the number of significantly differentially expressed genes
sum(resInvOffT0$padj < 0.1, na.rm =T)
sum(resInvOffT0$pvalue < 0.05, na.rm =T)

# Make a counts table that is scaled by the size factors
scaledcounts = counts(dds, normalize=T)
head(scaledcounts)

#building heat map data
head(scaledcounts)
genes4heatmap<-resInvOffT0[resInvOffT0$pvalue <0.05 & !is.na(resInvOffT0$pvalue),]
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
plotPCA(vstCounts, intgroup="Origin") +
  ggtitle("Deduped: Origin") +
  theme(legend.position = "bottom")
ggsave("Plots/dedup_origin_PCA.pdf")
  
plotPCA(vstCounts, intgroup="Geno") +
  ggtitle("Deduped: Genotype") +
  theme(legend.position = "bottom") 
ggsave("Plots/dedup_genotype_PCA.pdf")

plotPCA(vstCounts, intgroup="Temp") +
  ggtitle("Deduped: Temp") +
  theme(legend.position = "bottom")
ggsave("Plots/dedup_Temp_PCA.pdf")

plotPCA(vstCounts, intgroup="Timepoint") +
  ggtitle("Deduped: Timepoint") +
  theme(legend.position = "bottom")
ggsave("Plots/dedup_Timepoint_PCA.pdf")

NamesToSubset<-conds[conds$Timepoint =="TF", "Sample"]
plotPCA(vstCounts[,colnames(vstCounts)%in%NamesToSubset], intgroup="Origin")


####PCAs for specific results sets####
i="TF"
j="T0"
Group<-conds[conds$Timepoint==i | conds$Timepoint==j,"Sample"]
pdf(file=paste0("Plots/",i,"svs",j,"s_PCA_deduped.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()

i="TF"
Group<-conds[conds$Timepoint==i,"Sample"]
pdf(file=paste0("Plots/", i,"s_PCA_deduped.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()

i="T5"
Group<-conds[conds$Timepoint==i,"Sample"]
pdf(file=paste0("Plots/",i,"s_PCA_deduped.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Temp")
dev.off()

i="T0"
Group<-conds[conds$Timepoint==i,"Sample"]
pdf(file=paste0("Plots", i,"s_PCA_deduped.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin_Temp_Timepoint")
dev.off()

i="T0"
j="MoteOff"
Group<-conds[conds$Temp==i | conds$Origin==j,"Sample"]
pdf(file=paste0("Plots/", i,"svs",j,"s_PCA_deduped.pdf"),7, 7)
plotPCA(vstCounts[,colnames(vstCounts) %in% Group], intgroup="Origin")
dev.off()




#######################
######Jakes's Code#####
# vsd_Chlamy_N <- vst(dds_Chlamy_N, blind=FALSE)
#calculating the angle between two points in degrees

dev.off() 

# Libraries 
library(ggplot2)
library(tidyverse)
library(cowplot)

# Create PCA df #
pcaData_Chlamy_N <- plotPCA(vstCounts, intgroup=c("Temp","Timepoint", "Origin"), returnData=TRUE)
# save this and calc difference between sample and control
# what is the size of reaction of particular treatment group to our stress, is that different between our 2 groups?
  # within each TP, In vs Off, within each group, "For T1-T5 what's the reaction size of our samples in the different temp treatments" 
  # point for every individual and then box for each TP and Temperature 
  # try by hand for one TP/ group, control value for individual at that TP
  # make single vector in same order of every sample in order of geno and temp 
    # tidyverse trick to match samples by certain variable? 

#Add control values into the df and then do subtraction


PrinPCA <- princomp(scaledcounts, cor = T)
Loadings <- loadings(PrinPCA)
sdev <- PrinPCA$sdev
propvar <- PrinPCA$sdev^2/sum(PrinPCA$sdev^2)
cumsum(propvar)


T1_30_In <- pcaData_Chlamy_N %>%
  filter(group == "30:T1:MoteIn")

T1_30_Off <- pcaData_Chlamy_N %>%
  filter(group == "30:T1:MoteOff")

T2_30_In <- pcaData_Chlamy_N %>%
  filter(group == "30:T2:MoteIn")

T2_30_Off <- pcaData_Chlamy_N %>%
  filter(group == "30:T2:MoteOff")

T3_30_In <- pcaData_Chlamy_N %>%
  filter(group == "30:T3:MoteIn")

T3_30_Off <- pcaData_Chlamy_N %>%
  filter(group == "30:T3:MoteOff")

T4_30_In <- pcaData_Chlamy_N %>%
  filter(group == "30:T4:MoteIn")

T4_30_Off <- pcaData_Chlamy_N %>%
  filter(group == "30:T4:MoteOff")

T5_30_In <- pcaData_Chlamy_N %>%
  filter(group == "30:T5:MoteIn")

T5_30_Off <- pcaData_Chlamy_N %>%
  filter(group == "30:T5:MoteOff")

shift_data <- pcaData_Chlamy_N %>%
  filter(Timepoint != "T0" & Timepoint != "TF") %>%
  arrange(Timepoint, Origin, Temp) 

shift_data$ContPC1 <- c(rep(T1_30_In$PC1, 4), rep(T1_30_Off$PC1, 4),
                        rep(T2_30_In$PC1, 4), rep(T2_30_Off$PC1, 4),
                        rep(T3_30_In$PC1, 4), rep(T3_30_Off$PC1, 4),
                        rep(T4_30_In$PC1, 4), rep(T4_30_Off$PC1, 4),
                        rep(T5_30_In$PC1, 4), rep(T5_30_Off$PC1, 4))

shift_data$ContPC2 <- c(rep(T1_30_In$PC2, 4), rep(T1_30_Off$PC2, 4),
                        rep(T2_30_In$PC2, 4), rep(T2_30_Off$PC2, 4),
                        rep(T3_30_In$PC2, 4), rep(T3_30_Off$PC2, 4),
                        rep(T4_30_In$PC2, 4), rep(T4_30_Off$PC2, 4),
                        rep(T5_30_In$PC2, 4), rep(T5_30_Off$PC2, 4))

write.table(shift_data, file = "shift_data.txt", sep="\t", quote=F, row.names=TRUE)

percentVar_Chlamy_N <- round(100 * attr(pcaData_Chlamy_N, "percentVar"))

aggregate(PC1 ~ group, data=shift_data, length)

shift_data <- shift_data %>%
  mutate(PC1_shift = (PC1- ContPC1) * 0.30)

shift_data$PC1_shift <- abs(shift_data$PC1 - shift_data$ContPC1) * attr(pcaData_Chlamy_N, "percentVar")[1]

shift_data$PC2_shift <- abs(shift_data$PC2 - shift_data$ContPC2) * attr(pcaData_Chlamy_N, "percentVar")[2]


shift_data_shifted <- shift_data %>%
  dplyr::mutate(combo_shift = PC1_shift + PC2_shift) %>%
  filter(Temp != 30)

write.table(shift_data_shifted, file = "2024-07-09_shift_data_update.txt", sep="\t", quote=F, row.names=TRUE)


ggplot(shift_data_shifted, aes(x = Origin, y = combo_shift)) + #flip order of sites on axis so it goes from lowest to highest ed50
  geom_boxplot(aes(fill=Origin)) +
  scale_fill_manual(labels=c('Inshore', 'Offshore'), values=c('#8ADCE3', "#003A9D")) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white") +
  labs(y = "Transcriptomic Shift") +
  theme_bw() +
  facet_grid(Temp~Timepoint, scale = "free") +
  ylim(c(0,40)) +
  theme(legend.position = "none",
        # panel.spacing = unit(0, "lines"),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        # panel.grid.minor.y = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16))

ggplot(shift_data_shifted, aes(x = combo_shift)) +
  geom_histogram(binwidth = 20)

library(ggpubr)
library(car)

ggdensity(shift_data_shifted$combo_shift,
          main = "Density plot of transcriptomic shift",
          xlab = "Shift") 
ggqqplot(shift_data_shifted$combo_shift)

shapiro.test(shift_data_shifted$combo_shift)

bartlett.test(combo_shift ~ group, data = as.matrix(shift_data_shifted))


 




  



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

labs <- c("Inshore", "Offshore")
names(labs) <- c("MoteIn", "MoteOff")

ggplot(pcaData, aes(PC1, PC2,  shape=Origin)) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.85, aes(fill = Temp)) + #ylim(-30,30)+xlim(-30,30)+
  #geom_text(aes(label = Temp), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#00A5FF")) +
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+            
  facet_grid(Origin~Timepoint, labeller = labeller(Origin = labs)) +
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none",
        axis.text.y = element_text(angle = 0, size = 18),
        axis.text.x = element_text(angle = 0, size = 18),
        axis.title = element_text(size = 19, face = "italic"),
        strip.text = element_text(size = 19),
        strip.background = element_rect(color = "black", size = 1))

ggsave("Plots/PCA_TP-Facet_deduped.png", width = 17, height = 10, units = "in", dpi = "retina")

# Inshore vs Offshore by Temperature 
ggplot(pcaData, aes(PC1, PC2, fill = factor(Timepoint), shape=factor(Origin))) +
  # geom_spoke(aes(angle = pcaData_Chlamy_N$angle), arrow=arrow(length = unit(.3, "cm"), type="closed"),
  # radius =5.75, stat = "identity", position = "identity", inherit.aes = T, alpha=.7 )+
  geom_point(size=7,  alpha=.85) + #ylim(-30,30)+xlim(-30,30)+
  geom_text(aes(label = Timepoint), size =2)+
  scale_shape_manual(values=c(21,22), name ="Origin")+
  #scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF")) +
  xlab(paste0("PC1: ",percentVar_Chlamy_N[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_Chlamy_N[2],"% variance")) + 
  coord_fixed()+ggtitle("")+
  theme_light()+            
  facet_grid(Origin~Temp) +
  theme(strip.text = element_text(size = 16),
        axis.text.y = element_text(angle = 0, size = 14),
        axis.text.x = element_text(angle = 0, size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none")
ggsave("Plots/PCA_Temp-Facet_deduped.png", width = 16, height = 11, units = "in", dpi = "retina")

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
  scale_fill_manual(values = c("#00A5FF", "#00B81F",  "#BB9D00", "#F8766D", "#C77CFF", "#FF61C3")) +
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



#Here's code to subset the significant ones for individual contigs:
sigres_allContsdFilter<-res[res$padj<0.05 & !is.na(res$padj),]

#Merging annotations (metadata is a dataframe of the annotation table) with sigresults using the shared contignames:
write.table(cbind(Annotations[Annotations$ContigName%in%sigres_allContsdFilter[,1],], sigres_allContsdFilter, scaledcounts[row.names(scaledcounts)%in%row.names(sigres_allContsdFilter),]), file="12_300contvs400cont.txt", sep="\t", row.names=FALSE, quote=FALSE)

#################################################
########Trying GO Analysis from Sanbomics########
#################################################

# BiocManager::install("clusterProfiler")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("PCAtools")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(PCAtools)

sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05 & sigs$baseMean > 50,]

genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)

fit <- plot(barplot(GO_results, showCategory = 15))

png("out.png", res = 250, width = 1400, height = 1800)
print(fit)
dev.off()

fit

#################################################
#################################################
#################################################


#################################################
############### GOSeq analysis ##################
#################################################

Annotations<-read.delim("data/Sequencing/Porites_astreoides/Pilot_Fail_New_Concat_Compiled/P_astreoides_Kenkel2013_totalannotated.txt")
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

EnGOS<-GOGOGO$category[p.adjust(GOGOGO$over_represented_pvalue, method="fdr")<.05]


GOGOGO %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")

#################################################
#################################################
#################################################

#In_30_T1

scaledcounts






