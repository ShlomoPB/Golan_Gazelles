library(spAbundance)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the NB model (better fit)
fit <- readRDS("gazelle_msNMix_fit_nb_full.rds")
data_list <- readRDS("gazelle_msNMix_data_list.rds")

cat("Loaded NB model\n")

# Function to create effect size plot (posterior means with significance indicators)
plot_effect_sizes <- function(cov_name) {
  param_names <- dimnames(fit$beta.samples)[[2]]
  species_names <- c('male', 'female', 'juvenile')
  df <- data.frame()

  for (i in 1:3) {
    param_pattern <- paste0(cov_name, '-', species_names[i])
    cov_idx <- grep(param_pattern, param_names)
    if (length(cov_idx) == 0) next

    param_samples <- fit$beta.samples[, cov_idx]
    mean_val <- mean(param_samples)
    ci_low <- quantile(param_samples, 0.05)  # 90% CI lower
    ci_high <- quantile(param_samples, 0.95) # 90% CI upper

    # Determine significance
    significance <- ifelse(ci_low > 0, "positive",
                          ifelse(ci_high < 0, "negative", "uncertain"))

    df <- rbind(df, data.frame(species = species_names[i], cov = cov_name,
                               mean = mean_val, ci_low = ci_low, ci_high = ci_high,
                               significance = significance))
  }

  # Create plot
  p <- ggplot(df, aes(x = species, y = mean, color = significance)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', alpha = 0.5) +
    scale_color_manual(values = c("positive" = "blue", "negative" = "red", "uncertain" = "gray50"),
                       labels = c("positive" = "Significant +\n(90% CI > 0)",
                                "negative" = "Significant -\n(90% CI < 0)",
                                "uncertain" = "Uncertain\n(90% CI includes 0)")) +
    labs(title = paste('Effect Sizes for', cov_name, 'on Gazelle Abundance'),
         y = 'Log Abundance Change per SD Increase', x = 'Species',
         color = 'Effect Strength') +
    theme_minimal() +
    theme(text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)) +
    coord_flip()  # Horizontal bars for better readability

  return(p)
}

# Plot key covariates
covs <- c('Wolf.pack.abundance.z', 'null.jackal.abundance.z', 'Minefields.z', 'Proximity.to.Water.z')

for (cov in covs) {
  cat("Creating effect size plot for", cov, "\n")
  p <- plot_effect_sizes(cov)
  filename <- paste0('effect_sizes_', cov, '_nb.png')
  ggsave(filename, p, width = 8, height = 5, dpi = 150)
  cat("Saved", filename, "\n")
}

cat("Effect size plots completed.\n")