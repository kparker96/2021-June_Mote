library(readxl)
setwd("C:/Users/kpark/OneDrive/Documents/Barshis_Lab")


a <- read_xlsx(path = "reseq_rename.xlsx")

b <- read_xlsx(path = "reseq_filenames_test.xlsx") 

c <- right_join(a, b, by = "New_workbook_name")

d <- read_xlsx(path = "Mote_master_sample_name_sheet.xlsx")

library(tidyverse)

e <- left_join(c,d, by = "Novogene_R#") %>%
  select(`Date_RNAseq/Frozen`, Site, Species, Genotype, Timepoint, Treatment, Treatment_Temp_Setpoint, `Mote_Field/Frozen_SampleName`, Workbook__full_sample_label, Workbook_sample_token_label, Preferred_nomenclature, New_workbook_name = New_workbook_name.y, RNASeq_R1, RNASeq_R2, RNATubeNumber, `Novogene_R#`, `Reseq_Novogene_#`, RNASeq_R1_reseq, RNASeq_R2_reseq)

e$RNATubeNumber <- as.character(e$RNATubeNumber)

library(writexl)

write_xlsx(e,"reseq_metadata.xlsx")