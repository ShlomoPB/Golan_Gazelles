library(spAbundance)
library(tidyverse)
library(bayesplot)
library(tidybayes)

# Load the fitted lfMsNMix model
out <- readRDS("gazelle_lfMsNMix_fit_nb_null_wolf.rds")

# Extract species names
species_names <- c("Male", "Female", "Juvenile")

# ============================================================================
# 1. COMMUNITY-LEVEL EFFECT PLOT (All predictors with 95% CI)
# ============================================================================

com_summary <- summary(out, level = "community")
com_abund <- com_summary$abund.mean.table

# Create data frame for plotting
com_effects <- data.frame(
  predictor = rownames(com_abund),
  mean = com_abund[, "Mean"],
  lower = com_abund[, "2.5%"],
  upper = com_abund[, "97.5%"]
) %>%
  filter(predictor != "(Intercept)") %>%
  mutate(
    predictor = factor(predictor, levels = rev(predictor)),
    significant = !(lower < 0 & upper > 0)
  )

png("community_effects_lfMsNMix.png", width = 10, height = 6, res = 300)
ggplot(com_effects, aes(x = predictor, y = mean, color = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"), guide = "none") +
  labs(
    title = "Community-Level Abundance Effects (lfMsNMix Model)",
    subtitle = "Latent Factor Multispecies N-Mixture with Null.Wolf.Abundance",
    x = "Predictor (standardized)",
    y = "Effect Size (log scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text = element_text(size = 10)
  )
dev.off()

cat("Saved: community_effects_lfMsNMix.png\n")

# ============================================================================
# 2. SPECIES-LEVEL EFFECT PLOT (Comparing species responses to key predictors)
# ============================================================================

species_summary <- summary(out, level = "species")
species_abund <- species_summary$abund.table

# Parse species-specific effects
species_effects_list <- list()

for (pred in c("Null.Wolf.Abundance.z", "null.jackal.abundance.z", 
               "Built.and.Roads.z", "Nature.Reserves.z")) {
  for (i in seq_along(species_names)) {
    row_name <- paste0(pred, "-", tolower(species_names[i]))
    if (row_name %in% rownames(species_abund)) {
      species_effects_list[[length(species_effects_list) + 1]] <- data.frame(
        predictor = pred,
        species = species_names[i],
        mean = species_abund[row_name, "Mean"],
        lower = species_abund[row_name, "2.5%"],
        upper = species_abund[row_name, "97.5%"]
      )
    }
  }
}

species_effects <- do.call(rbind, species_effects_list) %>%
  mutate(
    species = factor(species, levels = species_names),
    predictor = gsub(".z$", "", predictor),
    significant = !(lower < 0 & upper > 0)
  )

png("species_effects_comparison_lfMsNMix.png", width = 12, height = 7, res = 300)
ggplot(species_effects, aes(x = species, y = mean, fill = species, color = significant)) +
  geom_col(position = "dodge", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(0.9), width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~predictor, scales = "free_y", nrow = 2) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "darkred"), guide = "none") +
  scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e", "Juvenile" = "#2ca02c")) +
  labs(
    title = "Species-Level Abundance Effects (lfMsNMix Model)",
    subtitle = "Comparing responses across male, female, and juvenile gazelles",
    x = "Species",
    y = "Effect Size (log scale)",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )
dev.off()

cat("Saved: species_effects_comparison_lfMsNMix.png\n")

# ============================================================================
# 3. WOLF ABUNDANCE EFFECT - FOREST PLOT (95% CI by species)
# ============================================================================

wolf_effects <- species_effects %>%
  filter(predictor == "Null.Wolf.Abundance") %>%
  arrange(mean)

png("wolf_effect_forest_plot_lfMsNMix.png", width = 9, height = 5, res = 300)
ggplot(wolf_effects, aes(x = reorder(species, mean), y = mean)) +
  geom_point(size = 4, color = "darkblue") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 height = 0.3, linewidth = 1, color = "darkblue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  coord_flip() +
  labs(
    title = "Null.Wolf.Abundance Effect on Gazelle Abundance",
    subtitle = "Forest plot with 95% credible intervals (lfMsNMix model)",
    x = "Species",
    y = "Effect Size (log scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray90")
  )
dev.off()

cat("Saved: wolf_effect_forest_plot_lfMsNMix.png\n")

# ============================================================================
# 4. CREDIBLE INTERVAL WIDTH COMPARISON (Convergence/Precision)
# ============================================================================

com_effects_ci_width <- data.frame(
  predictor = rownames(com_abund),
  ci_width = com_abund[, "97.5%"] - com_abund[, "2.5%"],
  rhat = com_abund[, "Rhat"],
  ess = com_abund[, "ESS"]
) %>%
  filter(predictor != "(Intercept)") %>%
  mutate(predictor = factor(predictor, levels = rev(predictor)))

png("convergence_diagnostics_lfMsNMix.png", width = 12, height = 6, res = 300)
par(mfrow = c(1, 2))

# CI width
plot(com_effects_ci_width$ci_width, 
     names.arg = com_effects_ci_width$predictor,
     main = "Community Effect - 95% CI Width",
     xlab = "Predictor",
     ylab = "CI Width",
     las = 2,
     col = ifelse(com_effects_ci_width$ci_width > 2, "coral", "steelblue"),
     ylim = c(0, max(com_effects_ci_width$ci_width) * 1.1))

# Rhat
plot(com_effects_ci_width$rhat,
     main = "Convergence Diagnostic (Rhat)",
     xlab = "Predictor",
     ylab = "Rhat",
     las = 2,
     col = ifelse(com_effects_ci_width$rhat < 1.01, "steelblue", "coral"),
     ylim = c(0.99, max(com_effects_ci_width$rhat, 1.02) * 1.05))
abline(h = 1.01, lty = 2, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()

cat("Saved: convergence_diagnostics_lfMsNMix.png\n")

# ============================================================================
# 5. DETECTION PROBABILITY EFFECTS (by species)
# ============================================================================

det_summary <- species_summary$det.table

det_effects_list <- list()

for (pred in c("log_days_operational.z", "cattle_sightings.z")) {
  for (i in seq_along(species_names)) {
    row_name <- paste0(pred, "-", tolower(species_names[i]))
    if (row_name %in% rownames(det_summary)) {
      det_effects_list[[length(det_effects_list) + 1]] <- data.frame(
        predictor = pred,
        species = species_names[i],
        mean = det_summary[row_name, "Mean"],
        lower = det_summary[row_name, "2.5%"],
        upper = det_summary[row_name, "97.5%"]
      )
    }
  }
}

det_effects <- do.call(rbind, det_effects_list) %>%
  mutate(
    species = factor(species, levels = species_names),
    predictor = gsub(".z$", "", predictor),
    predictor = factor(predictor, levels = c("log_days_operational", "cattle_sightings"))
  )

png("detection_effects_lfMsNMix.png", width = 10, height = 6, res = 300)
ggplot(det_effects, aes(x = species, y = mean, fill = species)) +
  geom_col(position = "dodge", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(0.9), width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~predictor, scales = "free_y") +
  scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e", "Juvenile" = "#2ca02c")) +
  labs(
    title = "Detection Probability Effects (lfMsNMix Model)",
    subtitle = "Species-specific responses to sampling effort and livestock sightings",
    x = "Species",
    y = "Effect Size (logit scale)",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )
dev.off()

cat("Saved: detection_effects_lfMsNMix.png\n")

# ============================================================================
# 6. INTERCEPT COMPARISON ACROSS SPECIES (baseline abundance)
# ============================================================================

intercept_data <- data.frame(
  species = species_names,
  mean = c(
    species_abund["(Intercept)-male", "Mean"],
    species_abund["(Intercept)-female", "Mean"],
    species_abund["(Intercept)-juvenile", "Mean"]
  ),
  lower = c(
    species_abund["(Intercept)-male", "2.5%"],
    species_abund["(Intercept)-female", "2.5%"],
    species_abund["(Intercept)-juvenile", "2.5%"]
  ),
  upper = c(
    species_abund["(Intercept)-male", "97.5%"],
    species_abund["(Intercept)-female", "97.5%"],
    species_abund["(Intercept)-juvenile", "97.5%"]
  )
) %>%
  mutate(species = factor(species, levels = species_names))

png("baseline_abundance_intercepts_lfMsNMix.png", width = 9, height = 6, res = 300)
ggplot(intercept_data, aes(x = species, y = mean, fill = species)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, linewidth = 0.8, color = "black") +
  scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e", "Juvenile" = "#2ca02c")) +
  labs(
    title = "Baseline Abundance by Species (Intercepts)",
    subtitle = "Average abundance on log scale when all predictors = 0",
    x = "Species",
    y = "Mean Abundance (log scale)",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )
dev.off()

cat("Saved: baseline_abundance_intercepts_lfMsNMix.png\n")

# ============================================================================
# Summary Statistics Export
# ============================================================================

cat("\n===========================================\n")
cat("lfMsNMix Model Visualization Summary\n")
cat("===========================================\n\n")

cat("Community-Level Effects:\n")
print(com_effects)

cat("\n\nWolf Effect by Species:\n")
print(wolf_effects)

cat("\n\nBaseline Abundance (log scale):\n")
print(intercept_data)

cat("\n\n✓ All visualizations saved successfully!\n")
