library(spAbundance)
library(readxl)
library(ggplot2)
library(dplyr)

# Load data
fit <- readRDS("gazelle_msNMix_fit.rds")
coords <- read_excel("Camera Coordinates.xlsx")
names(coords) <- c("site", "lat", "lon")
coords$site <- as.integer(coords$site)

# Abundance means: samples x species x sites
mu_mean <- apply(fit$mu.samples, c(2,3), mean)

# Create data frame for plotting
abund_df <- data.frame()
species_names <- c("male", "female", "juvenile")
for (sp in 1:3) {
  df_sp <- data.frame(site = 1:60, species = species_names[sp], abundance = mu_mean[sp, ])
  df_sp <- left_join(df_sp, coords, by = "site")
  abund_df <- rbind(abund_df, df_sp)
}

# Plot maps
p <- ggplot(abund_df, aes(x = lon, y = lat, color = abundance, size = abundance)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~species) +
  labs(title = "Predicted Gazelle Abundance", x = "Longitude", y = "Latitude") +
  theme_minimal()

ggsave("gazelle_abundance_map.png", p, width = 10, height = 8)

print("Abundance map saved.")