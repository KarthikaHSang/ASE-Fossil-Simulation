############################################################
## 02_analyze_results.R
##
## WHAT THIS SCRIPT DOES
## ----------------------------------------------------------
## Reads the simulation results produced by
## 01_simulate_fossils.R and:
##   - summarises ASE error by heterogeneity and overall
##     preservation level
##   - tests whether preservation level and heterogeneity
##     significantly affect ASE error (ANOVA)
##   - runs Tukey-style pairwise post-hoc comparisons
##   - produces all figures and summary tables
##
## This script is fast (seconds, not hours) -- it does NOT
## re-run any simulations.
##
## INPUT
## ----------------------------------------------------------
## simulation_results_revised.csv, as produced by
## 01_simulate_fossils.R (or your own copy of it, e.g.
## "FINAL_simulation_results_8000.csv" -- update the file name
## in SECTION 1 below to match whatever you have).
##
## OUTPUT
## ----------------------------------------------------------
## Several .png figures and printed summary tables (Table 1-4).
############################################################

## ---- Packages ------------------------------------------------
# install.packages("emmeans")  # run once if not already installed

library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)


############################################################
## SECTION 1: Load and prepare data
############################################################

results_file <- "simulation_results_revised.csv"
res <- read.csv(results_file)

## Consistent ordering for both factors, used throughout
category_levels <- c("None", "Low", "Medium", "High", "Very High", "No fossils")
preservation_levels <- c("None", "Low", "Medium", "High")

res$category <- factor(res$category, levels = category_levels)
res$OverallPreservation <- factor(res$OverallPreservation, levels = preservation_levels)

## Absolute difference between the two preservation rates --
## an alternative, continuous measure of heterogeneity
res$heterogeneity <- abs(res$psi_A - res$psi_B)

cat("Rows:", nrow(res), " | Columns:", ncol(res), "\n")
cat("Simulations per scenario:\n")
print(table(res$category, res$OverallPreservation))


############################################################
## SECTION 2: Table 1 -- experimental design summary
############################################################

table1 <- unique(res[, c("category", "psi_A", "psi_B", "OverallPreservation")])
table1$n <- 500
names(table1) <- c(
  "Heterogeneity", "Psi_A", "Psi_B", "Overall_Preservation", "Replicates"
)

print(table1)
write.csv(table1, "Table1_experimental_design.csv", row.names = FALSE)


############################################################
## SECTION 3: Table 2 -- descriptive statistics of ASE error
############################################################

table2 <- res %>%
  group_by(category, OverallPreservation) %>%
  summarise(
    n = n(),
    Mean_ASE_Error = round(mean(error), 4),
    Median_ASE_Error = round(median(error), 4),
    SD = round(sd(error), 4),
    .groups = "drop"
  ) %>%
  arrange(OverallPreservation, category)

print(table2)
write.csv(table2, "Table2_descriptive_stats.csv", row.names = FALSE)


############################################################
## SECTION 4: Figure 1 -- ASE error by preservation heterogeneity
##
## Excludes OverallPreservation == "None", since that facet
## only ever contains the single "No fossils" control and adds
## no comparative information.
############################################################

fig1_data <- subset(res, OverallPreservation != "None")

ggplot(fig1_data, aes(x = category, y = error, fill = category)) +
  geom_boxplot() +
  facet_wrap(~OverallPreservation) +
  labs(
    x = "Preservation heterogeneity",
    y = "ASE error",
    fill = "Preservation heterogeneity"
  ) +
  theme_classic()

ggsave(
  "Figure1_ASE_error_by_preservation_heterogeneity.png",
  width = 9, height = 8, dpi = 300
)


############################################################
## SECTION 5: Figure 2 -- interaction plot
############################################################

interaction_data <- res %>%
  filter(category != "No fossils") %>%
  group_by(category, OverallPreservation) %>%
  summarise(mean_error = mean(error), .groups = "drop")

interaction_data$OverallPreservation <- factor(
  interaction_data$OverallPreservation,
  levels = c("Low", "Medium", "High")
)

ggplot(
  interaction_data,
  aes(x = category, y = mean_error, group = OverallPreservation, linetype = OverallPreservation)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  labs(
    x = "Preservation heterogeneity",
    y = "Mean ASE error",
    linetype = "Overall preservation"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 11),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11)
  )

ggsave(
  "Figure2_interaction_heterogeneity_overall_preservation.png",
  width = 8, height = 6, dpi = 300
)


############################################################
## SECTION 6: Figure 3 -- tree shape vs ASE error
############################################################

ggplot(res, aes(x = tree_shape, y = error, colour = category)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~OverallPreservation) +
  labs(
    x = "Tree shape (Colless index)",
    y = "ASE error",
    colour = "Preservation heterogeneity"
  ) +
  theme_classic()

ggsave("Figure3_tree_shape_vs_ASE_error.png", width = 10, height = 6, dpi = 300)


############################################################
## SECTION 7: Figure 4 -- state transitions vs ASE error
############################################################

ggplot(res, aes(x = state_transitions, y = error, colour = category)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~OverallPreservation) +
  labs(
    x = "Number of state transitions",
    y = "ASE error",
    colour = "Preservation heterogeneity"
  ) +
  theme_classic()

ggsave("Figure4_state_transitions_vs_ASE_error.png", width = 10, height = 6, dpi = 300)


############################################################
## SECTION 8: Isolating the effect of overall preservation
##
## Restrict to simulations with NO heterogeneity (category ==
## "None", i.e. psi_A == psi_B), plus the no-fossil control
## (OverallPreservation == "None"). This isolates how much
## overall preservation intensity matters, independent of
## heterogeneity.
############################################################

preservation_none <- res %>%
  filter(category == "None" | OverallPreservation == "None")

preservation_none$total_preservation <- factor(
  preservation_none$OverallPreservation,
  levels = preservation_levels
)

summary_preservation <- preservation_none %>%
  group_by(total_preservation) %>%
  summarise(
    n = n(),
    mean_error = mean(error, na.rm = TRUE),
    sd_error = sd(error, na.rm = TRUE),
    se_error = sd_error / sqrt(n),
    CI_lower = mean_error - qt(0.975, df = n - 1) * se_error,
    CI_upper = mean_error + qt(0.975, df = n - 1) * se_error,
    .groups = "drop"
  )

## Magnitude of change relative to the no-fossil baseline
none_error <- summary_preservation$mean_error[
  summary_preservation$total_preservation == "None"
]

summary_preservation <- summary_preservation %>%
  mutate(
    absolute_reduction = none_error - mean_error,
    percentage_reduction = (none_error - mean_error) / none_error * 100
  )

print(summary_preservation)
write.csv(summary_preservation, "Table_preservation_effect.csv", row.names = FALSE)

ggplot(summary_preservation, aes(x = total_preservation, y = mean_error)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  labs(x = "Total preservation level", y = "Mean ASE error") +
  theme_classic()

ggsave("Figure5_total_preservation_vs_ASE_error.png", width = 7, height = 5, dpi = 300)

## One-way ANOVA: does overall preservation level affect error,
## when heterogeneity is held at zero?
preservation_anova <- aov(error ~ total_preservation, data = preservation_none)
cat("\n--- One-way ANOVA: preservation level (heterogeneity = None) ---\n")
print(summary(preservation_anova))


############################################################
## SECTION 9: Table 3 -- factorial ANOVA
##
## Tests the effects of heterogeneity, overall preservation,
## and their interaction, using only scenarios that actually
## have fossils (excludes the "No fossils"/"None" control,
## which has no heterogeneity level to compare).
############################################################

factorial_data <- res %>%
  filter(OverallPreservation != "None", category != "No fossils")

factorial_data$category <- factor(
  factorial_data$category,
  levels = c("None", "Low", "Medium", "High", "Very High")
)
factorial_data$OverallPreservation <- factor(
  factorial_data$OverallPreservation,
  levels = c("Low", "Medium", "High")
)

factorial_model <- lm(error ~ category * OverallPreservation, data = factorial_data)
factorial_anova <- anova(factorial_model)

cat("\n--- Factorial ANOVA: heterogeneity x overall preservation ---\n")
print(factorial_anova)

table3 <- data.frame(
  Effect = c(
    "Preservation heterogeneity",
    "Overall preservation",
    "Heterogeneity x overall preservation",
    "Residuals"
  ),
  df = c(
    factorial_anova["category", "Df"],
    factorial_anova["OverallPreservation", "Df"],
    factorial_anova["category:OverallPreservation", "Df"],
    factorial_anova["Residuals", "Df"]
  ),
  F_value = c(
    round(factorial_anova["category", "F value"], 3),
    round(factorial_anova["OverallPreservation", "F value"], 3),
    round(factorial_anova["category:OverallPreservation", "F value"], 3),
    NA
  ),
  P_value = c(
    factorial_anova["category", "Pr(>F)"],
    factorial_anova["OverallPreservation", "Pr(>F)"],
    factorial_anova["category:OverallPreservation", "Pr(>F)"],
    NA
  )
)

table3$P_value <- ifelse(
  is.na(table3$P_value), NA,
  ifelse(table3$P_value < 0.001, "<0.001", round(table3$P_value, 3))
)

print(table3)
write.csv(table3, "Table3_factorial_ANOVA.csv", row.names = FALSE)


############################################################
## SECTION 10: Table 4 -- Tukey-adjusted pairwise comparisons
##
## Uses emmeans to compare every pair of heterogeneity levels
## within each overall preservation level, with a Tukey
## adjustment for multiple comparisons.
############################################################

posthoc <- emmeans(
  factorial_model,
  pairwise ~ category | OverallPreservation,
  adjust = "tukey"
)

table4 <- as.data.frame(posthoc$contrasts)

table4$Significant <- ifelse(table4$p.value < 0.05, "Yes", "No")

table4$estimate <- round(table4$estimate, 4)
table4$SE <- round(table4$SE, 5)
table4$t.ratio <- round(table4$t.ratio, 3)
table4$p.value <- ifelse(
  table4$p.value < 0.001, "<0.001", sprintf("%.4f", table4$p.value)
)

names(table4) <- c(
  "Comparison", "Overall_Preservation", "Estimate",
  "SE", "df", "t_ratio", "P_value", "Significant"
)

print(table4)
write.csv(table4, "Table4_posthoc_comparisons.csv", row.names = FALSE)


############################################################
## SECTION 11: Figure 6 -- significance heatmap of pairwise comparisons
############################################################

tukey_data <- table4 %>%
  separate(Comparison, into = c("Het1", "Het2"), sep = " - ") %>%
  mutate(
    Het1 = factor(Het1, levels = category_levels),
    Het2 = factor(Het2, levels = category_levels),
    Overall_Preservation = factor(Overall_Preservation, levels = c("Low", "Medium", "High"))
  )

ggplot(tukey_data, aes(x = Het1, y = Het2, fill = Significant)) +
  geom_tile(color = "white") +
  facet_wrap(~Overall_Preservation) +
  scale_fill_manual(values = c("Yes" = "#2c7fb8", "No" = "#f0f0f0")) +
  labs(
    x = "Preservation heterogeneity",
    y = "Preservation heterogeneity",
    fill = "Significant\n(P < 0.05)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure6_posthoc_significance_heatmap.png", width = 9, height = 6, dpi = 300)

message("Analysis complete. Figures and tables saved to the working directory.")
