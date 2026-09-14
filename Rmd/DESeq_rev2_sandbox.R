####specific contrasts####
#### InvOff_34_T1 ####
date()
dds$Origin_Temp_Timepoint <- relevel(dds$Origin_Temp_Timepoint, ref = "MoteOff_34_T1")

ddsOff_34_T1 <- DESeq(dds)

dim(ddsOff_34_T1)
resultsNames(ddsOff_34_T1)

# Take the difference between Offshore 34°C T1 and Inshore 30°C T1, and subtract from it the differences between Offshore 30°C T1 and Inshore 30°C T1 and between Inshore 34°C T1 and Inshore 30°C T1.
# (Offshore 34°C T1 − Offshore 30°C T1) − (Inshore 34°C T1 − Inshore 30°C T1)


resOff_34_T1vIn_34_T1 <- results(ddsOff_34_T1, contrast = list(c("Origin_Temp_Timepoint_MoteOff_34_T1_vs_MoteIn_30_T1"),
                                                               c("Origin_Temp_Timepoint_MoteOff_30_T1_vs_MoteIn_30_T1",
                                                                 "Origin_Temp_Timepoint_MoteIn_34_T1_vs_MoteIn_30_T1")))
                                 
write.table(resOff_34_T1vIn_34_T1, file="data/Sequencing/Porites_astreoides/DESeq_rev2/resOff_34_T1vIn_34_T1.txt", quote=F, row.names=T, sep="\t")
head(resOff_34_T1vIn_34_T1)
summary(resOff_34_T1vIn_34_T1)
date()