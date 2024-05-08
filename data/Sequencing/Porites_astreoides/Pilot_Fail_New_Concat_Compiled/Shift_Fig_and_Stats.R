library(tidyverse)
library(ggplot2)
library(car)

load("~/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing/Porites_astreoides/Pilot_Fail_New_Concat_Compiled/DESeq_dds.RData")

# TRANSCRIPTOMIC SHIFT DATA # 
shift_data_shifted <- read.delim("shift_data_shifted.txt") %>%
  dplyr::mutate(combo_shift = PC1_shift + PC2_shift) %>%
  filter(Temp != 30)

ggplot(shift_data_shifted, aes(x = Origin, y = combo_shift)) + #flip order of sites on axis so it goes from lowest to highest ed50
  geom_boxplot(aes(fill=Origin)) +
  scale_fill_manual(labels=c('Inshore', 'Offshore'), values=c('#8ADCE3', "#003A9D")) +
  stat_summary(fun = mean, geom = "point", shape = 23, fill = "white") +
  labs(y = "Transcriptomic Shift") +
  theme_bw() +
  facet_grid(Temp~Timepoint, scale = "free") +
  ylim(c(0,40)) +
  theme(legend.position = "none",
        panel.spacing = unit(0, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16))

ggplot(shift_data_shifted, aes(x = combo_shift)) +
  geom_histogram(binwidth = 20)

# STATS #

# Load your data (assuming your data frame is named "shift_data_shifted")

# 1. Check for independence of observations (nothing to code, already addressed).

# 2. Normality assumption (Shapiro-Wilk test)
shapiro_test <- shapiro.test(shift_data_shifted$combo_shift)

print("Shapiro-Wilk Normality Test:")

print(shapiro_test)

if (shapiro_test$p.value < 0.05) {
  print("The data may not be normally distributed.")
} else {
  print("The data appears to be normally distributed.")
}

# 3. Homogeneity of variances (Levene's test)
levene_test <- leveneTest(combo_shift ~ group, data = shift_data_shifted)

print("Levene's Test for Homogeneity of Variance:")

print(levene_test)

if (levene_test$p.value < 0.05) {
  print("There may be an issue with homogeneity of variances.")
} else {
  print("The variances appear to be homogeneous across groups.")
}

# 4. Independence of groups (nothing to code, already addressed).

# 5. Continuous dependent variable (nothing to code, you have a continuous variable).

# Residuals vs. Fitted plot (optional but useful)
model <- aov(combo_shift ~ group, data = shift_data_shifted)
plot(model)

# You can also create Q-Q plots to visualize normality.
qqnorm(resid(model))
qqline(resid(model))

aov(combo_shift ~ group,
    data = shift_data_shifted
)

oneway.test(combo_shift ~ group,
            data = shift_data_shifted,
            var.equal = TRUE # assuming equal variances
)

# Perform Tukey's HSD post hoc test
tukey_result <- TukeyHSD(model) 



