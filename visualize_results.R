library(spAbundance)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load fitted model and data
fit <- readRDS("gazelle_msNMix_fit.rds")
data_list <- readRDS("gazelle_msNMix_data_list.rds")

# Extract posterior samples
beta.samples <- fit$beta.samples  # Community abundance coefficients
alpha.samples <- fit$alpha.samples  # Community detection coefficients
beta.samples.sp <- fit$beta.samples  # Species-specific (wait, in msNMix, beta.samples is community, beta.samples.sp is species)
# Actually, for msNMix, beta.comm.samples, alpha.comm.samples, beta.samples (species), alpha.samples (species)

# From help: beta.comm.samples, alpha.comm.samples, beta.samples (array species x params x samples), etc.

# For simplicity, use summary for plots.

# Function to plot covariate effects
plot_cov_effect <- function(cov_name, species = NULL, level = 'community') {
  if (level == 'community') {
    summ <- summary(fit, level = 'community')
    abund <- summ$abund.means
    row <- grep(cov_name, rownames(abund))
    if (length(row) == 0) return(NULL)
    mean_val <- abund[row, 'Mean']
    sd_val <- abund[row, 'SD']
    ci_low <- abund[row, '2.5%']
    ci_high <- abund[row, '97.5%']
    df <- data.frame(cov = cov_name, mean = mean_val, sd = sd_val, ci_low = ci_low, ci_high = ci_high)
  } else {
    summ <- summary(fit, level = 'species')
    abund <- summ$abund.means
    rows <- grep(cov_name, rownames(abund))
    if (length(rows) == 0) return(NULL)
    species_names <- c('male', 'female', 'juvenile')
    df <- data.frame()
    for (i in 1:3) {
      row <- rows[i]
      df <- rbind(df, data.frame(species = species_names[i], cov = cov_name,
                                 mean = abund[row, 'Mean'], sd = abund[row, 'SD'],
                                 ci_low = abund[row, '2.5%'], ci_high = abund[row, '97.5%']))
    }
  }
  p <- ggplot(df, aes(x = cov, y = mean, ymin = ci_low, ymax = ci_high)) +
    geom_point() + geom_errorbar(width = 0.2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(title = paste('Effect of', cov_name, 'on Abundance'), y = 'Log Abundance Change') +
    theme_minimal()
  if (!is.null(species)) p <- p + facet_wrap(~species)
  return(p)
}

# Plot key covariates
covs <- c('Wolf.pack.abundance.z', 'null.jackal.abundance.z', 'Minefields.z', 'Proximity.to.Water.z')

for (cov in covs) {
  p <- plot_cov_effect(cov, level = 'species')
  ggsave(paste0('effect_', cov, '.png'), p, width = 8, height = 6)
}

# Abundance map if coordinates available
# Assume coords <- read.csv('Camera Coordinates.xlsx') or something
# For now, skip

print("Visualization script run. Check PNG files.")