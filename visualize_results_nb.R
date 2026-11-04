library(spAbundance)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the NB model (better fit)
fit <- readRDS("gazelle_msNMix_fit_nb_full.rds")
data_list <- readRDS("gazelle_msNMix_data_list.rds")

cat("Loaded NB model\n")
cat("beta.samples dimensions:", dim(fit$beta.samples), "\n")
cat("Parameter names:", dimnames(fit$beta.samples)[[2]], "\n")

# Function to plot covariate effects with 90% CI
plot_cov_effect <- function(cov_name, level = 'species', ci_level = 0.90) {
  if (level == 'community') {
    summ <- summary(fit, level = 'community')
    abund <- summ$abund.means
    row <- grep(cov_name, rownames(abund))
    if (length(row) == 0) {
      cat("Covariate", cov_name, "not found in community summary\n")
      return(NULL)
    }
    mean_val <- abund[row, 'Mean']
    sd_val <- abund[row, 'SD']
    ci_low <- abund[row, '2.5%']
    ci_high <- abund[row, '97.5%']
    df <- data.frame(cov = cov, mean = mean_val, sd = sd_val, ci_low = ci_low, ci_high = ci_high)
    p <- ggplot(df, aes(x = cov, y = mean, ymin = ci_low, ymax = ci_high)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
      labs(title = paste('Community Effect of', cov_name, 'on Abundance'),
           y = 'Log Abundance Change', x = '') +
      theme_minimal() +
      theme(text = element_text(size = 12))
  } else {
    # Find the covariate indices for each species
    param_names <- dimnames(fit$beta.samples)[[2]]
    species_names <- c('male', 'female', 'juvenile')
    df <- data.frame()
    for (i in 1:3) {
      param_pattern <- paste0(cov_name, '-', species_names[i])
      cov_idx <- grep(param_pattern, param_names)
      if (length(cov_idx) == 0) {
        cat("Parameter", param_pattern, "not found\n")
        next
      }
      # beta.samples is [samples, parameters]
      param_samples <- fit$beta.samples[, cov_idx]
      mean_val <- mean(param_samples)
      sd_val <- sd(param_samples)
      alpha <- (1 - ci_level) / 2
      ci_low <- quantile(param_samples, alpha)
      ci_high <- quantile(param_samples, 1 - alpha)
      
      df <- rbind(df, data.frame(species = species_names[i], cov = cov_name,
                                 mean = mean_val, sd = sd_val,
                                 ci_low = ci_low, ci_high = ci_high))
    }
    cat("Data frame for", cov_name, ":\n")
    print(df)
    p <- ggplot(df, aes(x = species, y = mean, ymin = ci_low, ymax = ci_high, color = species)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
      labs(title = paste('Species-Specific Effect of', cov_name, 'on Abundance'),
           y = 'Log Abundance Change', x = 'Species') +
      scale_color_manual(values = c('male' = 'blue', 'female' = 'red', 'juvenile' = 'green')) +
      theme_minimal() +
      theme(text = element_text(size = 12), legend.position = 'none')
  }
  return(p)
}

# Plot key covariates
covs <- c('Wolf.pack.abundance.z', 'null.jackal.abundance.z', 'Minefields.z', 'Proximity.to.Water.z')

for (cov in covs) {
  cat("Creating plot for", cov, "\n")
  p <- plot_cov_effect(cov, level = 'species', ci_level = 0.90)
  if (!is.null(p)) {
    filename <- paste0('effect_', cov, '_nb_90ci.png')
    cat("Attempting to save to", filename, "\n")
    ggsave(filename, p, width = 8, height = 6, dpi = 150)
    cat("Saved", filename, "\n")
  }
}

# Also create community-level plots
summ <- summary(fit, level = 'community')
abund <- summ$abund.means

for (cov in covs) {
  cat("Creating community plot for", cov, "\n")
  row <- grep(cov, rownames(abund))
  if (length(row) > 0) {
    mean_val <- abund[row, 'Mean']
    sd_val <- abund[row, 'SD']
    ci_low <- abund[row, '2.5%']
    ci_high <- abund[row, '97.5%']
    df <- data.frame(cov = cov, mean = mean_val, sd = sd_val, ci_low = ci_low, ci_high = ci_high)
    p <- ggplot(df, aes(x = cov, y = mean, ymin = ci_low, ymax = ci_high)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
      labs(title = paste('Community Effect of', cov, 'on Abundance'),
           y = 'Log Abundance Change', x = '') +
      theme_minimal() +
      theme(text = element_text(size = 12))
    filename <- paste0('community_effect_', cov, '_nb.png')
    ggsave(filename, p, width = 6, height = 4, dpi = 150)
    cat("Saved", filename, "\n")
  } else {
    cat("Covariate", cov_comm, "not found in community summary\n")
  }
}

cat("Visualization script completed. Check PNG files.\n")