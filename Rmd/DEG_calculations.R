library(dplyr)

# ---------------------------------------------------------
# 1. Set folder
# ---------------------------------------------------------

folder <- "data/Sequencing/Porites_astreoides/DESeq_rev2/Reviewer-Code-Run/39"

# Timepoints
timepoints <- paste0("T", 1:5)

# Empty list to store results
all_results <- list()


# ---------------------------------------------------------
# 2. Read and classify each timepoint
# ---------------------------------------------------------

for (tp in timepoints) {
  
  # Read interaction results
  interaction <- read.delim(
    file.path(
      folder,
      paste0("rev2_res_interaction-test_39_", tp, ".txt")
    ),
    row.names = 1
  )
  
  # Read Offshore heat-response results
  offshore <- read.delim(
    file.path(
      folder,
      paste0("rev2_res_Offshore_39-30_", tp, ".txt")
    ),
    row.names = 1
  )
  
  # Read Inshore heat-response results
  inshore <- read.delim(
    file.path(
      folder,
      paste0("rev2_res_Inshore_39-30_", tp, ".txt")
    ),
    row.names = 1
  )
  
  
  # -------------------------------------------------------
  # 3. Keep only significant interaction genes
  # -------------------------------------------------------
  
  interaction_sig <- interaction %>%
    filter(!is.na(padj), padj < 0.1)
  
  
  # -------------------------------------------------------
  # 4. Combine the three results
  # -------------------------------------------------------
  
  combined <- data.frame(
    Gene = rownames(interaction_sig),
    Interaction_LFC = interaction_sig$log2FoldChange,
    Interaction_padj = interaction_sig$padj,
    Offshore_LFC = offshore[rownames(interaction_sig), "log2FoldChange"],
    Inshore_LFC = inshore[rownames(interaction_sig), "log2FoldChange"],
    Timepoint = tp,
    row.names = NULL
  )
  
  
  # -------------------------------------------------------
  # 5. Classify genes
  # -------------------------------------------------------
  
  combined <- combined %>%
    mutate(
      Category = case_when(
        
        # Offshore response only
        Offshore_LFC != 0 & Inshore_LFC == 0 ~
          "Offshore only",
        
        # Inshore response only
        Offshore_LFC == 0 & Inshore_LFC != 0 ~
          "Inshore only",
        
        # Opposite directions
        sign(Offshore_LFC) != sign(Inshore_LFC) ~
          "Opposite direction",
        
        # Same direction, stronger Offshore
        abs(Offshore_LFC) > abs(Inshore_LFC) ~
          "Same direction, stronger Offshore",
        
        # Same direction, stronger Inshore
        abs(Inshore_LFC) > abs(Offshore_LFC) ~
          "Same direction, stronger Inshore",
        
        # Anything unexpected
        TRUE ~ "Check"
      )
    )
  
  
  # Store this timepoint
  all_results[[tp]] <- combined
}


# ---------------------------------------------------------
# 6. Combine all timepoints
# ---------------------------------------------------------

classification_results <- bind_rows(all_results)

# ---------------------------------------------------------
# 7. Count categories
# ---------------------------------------------------------

category_counts <- classification_results %>%
  count(Timepoint, Category)

# View counts
category_counts

category_counts_same_direction <- category_counts %>%
  filter(Category == "Same direction, stronger Inshore" | Category == "Same direction, stronger Offshore")

library(ggplot2)

ggplot(category_counts_same_direction, aes(x = Timepoint, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=c('#8ADCE3', "#003A9D")) +
  ylab("# of DEGs from Interaction Test") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle("39C Treatment")

category_counts_unique <- category_counts %>%
  filter(Category == "Inshore only" | Category == "Offshore only")

ggplot(category_counts_unique, aes(x = Timepoint, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=c('#8ADCE3', "#003A9D")) +
  theme(legend.position = "bottom")

class_opposite_driections <- classification_results %>%
  filter(Category == "Opposite direction") 

class_opposite_long <- class_opposite_driections %>%
  select(Gene, Timepoint, Offshore_LFC, Inshore_LFC) %>%
  pivot_longer(
    cols = c(Offshore_LFC, Inshore_LFC),
    names_to = "Population",
    values_to = "LFC"
  ) %>%
  mutate(
    Population = recode(
      Population,
      Offshore_LFC = "Offshore",
      Inshore_LFC = "Inshore"
    )
  )

ggplot(class_opposite_long, aes(x = Gene, y = LFC, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values=c('#8ADCE3', "#003A9D")) +
  facet_wrap(~Timepoint, scales = "free_x") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank()) 
  
  
# ---------------------------------------------------------
# 8. View any remaining "Check" genes
# ---------------------------------------------------------

classification_results %>%
  filter(Category == "Check") %>%
  select(
    Gene,
    Timepoint,
    Offshore_LFC,
    Inshore_LFC
  )


# ---------------------------------------------------------
# 9. Save full classification table
# ---------------------------------------------------------

write.table(
  classification_results,
  file = file.path(
    folder,
    "rev2_39C_interaction_classification_all_timepoints.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)


# ---------------------------------------------------------
# 10. Save category counts
# ---------------------------------------------------------

write.table(
  category_counts,
  file = file.path(
    folder,
    "rev2_39C_interaction_category_counts.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)