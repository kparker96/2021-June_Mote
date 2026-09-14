load("Rmd/dds.Rdata")


vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)

genes_to_check <- c(
  "isotig04273_Past",
  "isotig01060_Past",
  "isotig01034_Past",
  "isotig07479_Past",
  "isotig07883_Past"
)

vst_check <- vst_counts[genes_to_check, ]

colData(dds)[, "Origin_Temp_Timepoint"]

group <- as.character(colData(dds)$Origin_Temp_Timepoint)

group_T3 <- group %in% c("MoteOff_34_T3", "MoteIn_34_T3")

vst_T3 <- vst_counts[genes_to_check, group_T3]

counts_check <- countsTable[genes_to_check, 
                            colnames(vst_T3)]

counts_check

library(tidyverse)

plot_df <- as.data.frame(vst_T3) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "VST"
  ) %>%
  mutate(
    Origin = ifelse(
      str_detect(Sample, "MoteOff"),
      "Offshore",
      "Inshore"
    )
  )

head(plot_df)

ggplot(plot_df, aes(x = Origin, y = VST, color = Origin)) +
  geom_jitter(width = 0.1, size = 2.5, alpha = 0.8) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 18,
    color = "black",
  ) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(
    x = NULL,
    y = "VST expression"
  ) +
  ggtitle("Offshore Downregulated genes at 34 T3") +
  theme_classic()