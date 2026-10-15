####load libraries####
library(tidyverse)
library(readxl)
library(writexl)

###additional DESeq Libraries###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
#BiocManager::install("geneLenDataBase")
#BiocManager::install("goseq")
#BiocManager::install("gplots")

library(DESeq2) 
library(gplots)

# setwd("/Users/kpark/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing/")  # The drag and drop from finder works in R, too.

####read in your data to make counts table####
countsTable_raw <- read_delim('data/Sequencing/Porites_astreoides/Counts/STAR_Salmon_Past_Smic_NumReads_matrix_Past.txt')
dim(countsTable_raw)

###read in sample names to change novogene_# to coral id number)
sample_names <- read_excel('data/Sequencing/Porites_astreoides/Renamer/Sample_Labels_Novogene_vs_RNA#.xlsx') %>%
  dplyr::select(sample = `Novogene_R#`, coral_ID = 'Mote_Field/Frozen_SampleName')

####rename columns to match sample names####
countsTable_clean <- countsTable_raw %>%
  gather(key = "sample", value = 'reads', -ContigName)

countsTable_newname <- right_join(countsTable_clean, sample_names, by = "sample") %>%
  dplyr::select(ContigName, sample = coral_ID, reads) %>%
  filter(reads >= 0)

countsTable <- pivot_wider(countsTable_newname, names_from = sample, values_from = reads) %>%
  column_to_rownames(var="ContigName")

colSums(countsTable)
colSums(countsTable)[order(colSums(countsTable), decreasing=T)]

write.table(countsTable, "data/Sequencing/Porites_astreoides/DESeq_rev2/Dedup_Mote_Past_Counts_full_cleannames.txt")

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

####check replications####
table(conds$Origin)
table(conds$Geno)
table(conds$Temp)
table(conds$Timepoint)
table(conds$Origin_Temp_Timepoint)



#### how many genes have mean count >=3 ####
#1 → row-wise (genes)
means = apply(countsTable, 1, mean)
table(means>=3)

# Filter low expressed contigs 
# Keep only contigs where the average count across all 216 samples is ≥ 3 
countsTable<-countsTable[means>=3,] 
countsTable<-round(countsTable, digits=0)
dim(countsTable)

total_reads <- sum(countsTable, na.rm = TRUE)
n_contigs <- nrow(countsTable)

library_reads <- colSums(countsTable, na.rm = TRUE)
min_reads <- min(library_reads)
max_reads <- max(library_reads)
mean_reads <- mean(library_reads)

#### count data sets####
load("Rmd/dds.Rdata")

#### Reviewer interaction contrast ####

# Check that the expected factor and coefficients exist
levels(colData(dds)$Origin_Temp_Timepoint)
resultsNames(dds)

# Define reviewer function
group_contrast <- function(dds, weights, factor = "Origin_Temp_Timepoint") {
  stopifnot(abs(sum(weights)) < 1e-8)
  
  rn <- resultsNames(dds)
  ref <- levels(colData(dds)[[factor]])[1]
  
  v <- setNames(numeric(length(rn)), rn)
  
  for (g in names(weights)) {
    if (g == ref) next
    
    coef_name <- paste0(factor, "_", g, "_vs_", ref)
    
    if (!coef_name %in% rn) {
      stop("Coefficient not found: ", coef_name)
    }
    
    v[coef_name] <- weights[[g]]
  }
  
  v
}

#### Interaction: (Offshore 37 - Offshore 30) - (Inshore 37 - Inshore 30), T1 ####

w <- c(
  MoteOff_37_T1 = 1,
  MoteOff_30_T1 = -1,
  MoteIn_37_T1 = -1,
  MoteIn_30_T1 = 1
)

res_int_T1 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

# Look at results
head(res_int_T1)
summary(res_int_T1)

# Save results
write.table(
  res_int_T1,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_interaction-test_37_T1.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)

#### Interaction T2 ####

w <- c(
  MoteOff_37_T2 = 1,
  MoteOff_30_T2 = -1,
  MoteIn_37_T2 = -1,
  MoteIn_30_T2 = 1
)

res_int_T2 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_int_T2,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_interaction-test_37_T2.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### Interaction T3 ####

w <- c(
  MoteOff_37_T3 = 1,
  MoteOff_30_T3 = -1,
  MoteIn_37_T3 = -1,
  MoteIn_30_T3 = 1
)

res_int_T3 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_int_T3,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_interaction-test_37_T3.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### Interaction T4 ####

w <- c(
  MoteOff_37_T4 = 1,
  MoteOff_30_T4 = -1,
  MoteIn_37_T4 = -1,
  MoteIn_30_T4 = 1
)

res_int_T4 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_int_T4,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_interaction-test_37_T4.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### Interaction T5 ####

w <- c(
  MoteOff_37_T5 = 1,
  MoteOff_30_T5 = -1,
  MoteIn_37_T5 = -1,
  MoteIn_30_T5 = 1
)

res_int_T5 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_int_T5,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_interaction-test_37_T5.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)

#### Population-specific heat contrasts ####
#### Offshore 37 - 30 and Inshore 37 - 30 ####


#### T1 ####

# Offshore 37 - 30
w <- c(
  MoteOff_37_T1 = 1,
  MoteOff_30_T1 = -1
)

res_off_T1 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_off_T1,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Offshore_37-30_T1.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


# Inshore 37 - 30
w <- c(
  MoteIn_37_T1 = 1,
  MoteIn_30_T1 = -1
)

res_in_T1 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_in_T1,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Inshore_37-30_T1.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### T2 ####

# Offshore 37 - 30
w <- c(
  MoteOff_37_T2 = 1,
  MoteOff_30_T2 = -1
)

res_off_T2 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_off_T2,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Offshore_37-30_T2.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


# Inshore 37 - 30
w <- c(
  MoteIn_37_T2 = 1,
  MoteIn_30_T2 = -1
)

res_in_T2 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_in_T2,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Inshore_37-30_T2.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### T3 ####

# Offshore 37 - 30
w <- c(
  MoteOff_37_T3 = 1,
  MoteOff_30_T3 = -1
)

res_off_T3 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_off_T3,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Offshore_37-30_T3.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


# Inshore 37 - 30
w <- c(
  MoteIn_37_T3 = 1,
  MoteIn_30_T3 = -1
)

res_in_T3 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_in_T3,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Inshore_37-30_T3.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### T4 ####

# Offshore 37 - 30
w <- c(
  MoteOff_37_T4 = 1,
  MoteOff_30_T4 = -1
)

res_off_T4 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_off_T4,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Offshore_37-30_T4.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


# Inshore 37 - 30
w <- c(
  MoteIn_37_T4 = 1,
  MoteIn_30_T4 = -1
)

res_in_T4 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_in_T4,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Inshore_37-30_T4.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


#### T5 ####

# Offshore 37 - 30
w <- c(
  MoteOff_37_T5 = 1,
  MoteOff_30_T5 = -1
)

res_off_T5 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_off_T5,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Offshore_37-30_T5.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)


# Inshore 37 - 30
w <- c(
  MoteIn_37_T5 = 1,
  MoteIn_30_T5 = -1
)

res_in_T5 <- results(
  dds,
  contrast = group_contrast(dds, w),
  alpha = 0.1
)

write.table(
  res_in_T5,
  file = "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/rev2_res_Inshore_37-30_T5.txt",
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)
