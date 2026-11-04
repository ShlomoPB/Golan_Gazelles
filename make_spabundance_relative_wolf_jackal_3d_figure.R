#!/usr/bin/env Rscript

# Static 3D companion to the publication wolf-jackal heatmap. This script uses
# saved posterior draws from the final spAbundance models and does not refit.

fit_file <- "publication_wolf_jackal_interaction_nmix_spabundance_fits.rds"
reference_grid_file <-
  "publication_spabundance_wolf_jackal_relative_prediction_grid.csv"
output_prefix <- "publication_spabundance_wolf_jackal_relative_3d"

class_labels <- c(fawn = "Juvenile", female = "Female", male = "Male")
predator_levels <- seq(-1, 1, length.out = 41)
reference_levels <- c(-1, 0, 1)

low_color <- "#C85C4A"
mid_color <- "#F6F3EA"
high_color <- "#204D69"
surface_palette <- grDevices::colorRampPalette(
  c(low_color, mid_color, high_color),
  space = "Lab"
)(256)

required_beta_names <- c(
  "Null.Wolf.Abundance.z",
  "null.jackal.abundance.z",
  "Null.Wolf.Abundance.z:null.jackal.abundance.z"
)

summarize_surface <- function(fit, grid, chunk_size = 100L) {
  draws <- as.matrix(fit$beta.samples)
  missing_beta <- setdiff(required_beta_names, colnames(draws))
  if (length(missing_beta) > 0L) {
    stop("Missing abundance coefficient(s): ", paste(missing_beta, collapse = ", "))
  }

  beta <- draws[, required_beta_names, drop = FALSE]
  design <- cbind(
    grid$wolf_z,
    grid$jackal_z,
    grid$wolf_z * grid$jackal_z
  )

  estimates <- numeric(nrow(grid))
  conf_low <- numeric(nrow(grid))
  conf_high <- numeric(nrow(grid))
  chunks <- split(
    seq_len(nrow(grid)),
    ceiling(seq_len(nrow(grid)) / chunk_size)
  )

  for (index in chunks) {
    relative_draws <- exp(beta %*% t(design[index, , drop = FALSE]))
    estimates[index] <- apply(relative_draws, 2, stats::median)
    conf_low[index] <- apply(
      relative_draws,
      2,
      stats::quantile,
      probs = 0.025,
      names = FALSE
    )
    conf_high[index] <- apply(
      relative_draws,
      2,
      stats::quantile,
      probs = 0.975,
      names = FALSE
    )
  }

  data.frame(
    grid,
    estimate = estimates,
    conf_low = conf_low,
    conf_high = conf_high,
    check.names = FALSE
  )
}

fits <- readRDS(fit_file)
if (!all(names(class_labels) %in% names(fits))) {
  stop("The saved fit object does not contain all three demographic classes.")
}

surface_grid <- expand.grid(
  jackal_z = predator_levels,
  wolf_z = predator_levels,
  KEEP.OUT.ATTRS = FALSE
)
reference_grid <- expand.grid(
  jackal_z = reference_levels,
  wolf_z = reference_levels,
  KEEP.OUT.ATTRS = FALSE
)

surface_results <- do.call(
  rbind,
  lapply(names(class_labels), function(class_name) {
    result <- summarize_surface(fits[[class_name]], surface_grid)
    result$class <- unname(class_labels[class_name])
    result
  })
)
row.names(surface_results) <- NULL

reference_results <- do.call(
  rbind,
  lapply(names(class_labels), function(class_name) {
    result <- summarize_surface(fits[[class_name]], reference_grid)
    result$class <- unname(class_labels[class_name])
    result
  })
)
row.names(reference_results) <- NULL

# Guard against drawing from a stale or different fit object.
saved_reference <- utils::read.csv(reference_grid_file, check.names = FALSE)
check_columns <- c("class", "wolf_z", "jackal_z", "estimate")
if (!all(check_columns %in% names(saved_reference))) {
  stop("The saved reference grid is missing required verification columns.")
}
verification <- merge(
  reference_results[, check_columns],
  saved_reference[, check_columns],
  by = c("class", "wolf_z", "jackal_z"),
  suffixes = c("_recomputed", "_saved"),
  sort = TRUE
)
if (
  nrow(verification) != 27L ||
    max(abs(verification$estimate_recomputed - verification$estimate_saved)) > 1e-10
) {
  stop("Recomputed reference predictions do not match the published heatmap grid.")
}

utils::write.csv(
  surface_results,
  paste0(output_prefix, "_surface.csv"),
  row.names = FALSE,
  na = ""
)

log_range <- range(log(surface_results$estimate))
z_upper <- ceiling(max(surface_results$estimate) * 4) / 4
z_limits <- c(0, z_upper)

rescale_log_mid <- function(values) {
  log_values <- log(values)
  scaled <- ifelse(
    log_values <= 0,
    0.5 * (log_values - log_range[1]) / (0 - log_range[1]),
    0.5 + 0.5 * log_values / log_range[2]
  )
  pmin(1, pmax(0, scaled))
}

abundance_colors <- function(values) {
  color_index <- 1L + floor(rescale_log_mid(values) * 255)
  surface_palette[color_index]
}

draw_color_key <- function() {
  graphics::par(mar = c(4.5, 0.2, 4.4, 3.6))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = log_range)

  key_edges <- seq(log_range[1], log_range[2], length.out = 257)
  key_midpoints <- (key_edges[-1] + key_edges[-length(key_edges)]) / 2
  graphics::rect(
    xleft = 0.08,
    ybottom = key_edges[-length(key_edges)],
    xright = 0.42,
    ytop = key_edges[-1],
    col = abundance_colors(exp(key_midpoints)),
    border = NA
  )

  key_breaks <- c(0.25, 0.5, 1, 2)
  key_breaks <- key_breaks[
    log(key_breaks) >= log_range[1] & log(key_breaks) <= log_range[2]
  ]
  graphics::axis(
    side = 4,
    at = log(key_breaks),
    labels = format(key_breaks, trim = TRUE),
    las = 1,
    tck = -0.18,
    mgp = c(2, 0.45, 0),
    cex.axis = 0.82
  )
  graphics::mtext(
    "Gazelle\nrelative\nabundance",
    side = 3,
    line = 0.2,
    adj = 0.5,
    cex = 0.76,
    font = 2
  )
}

draw_surface_panel <- function(class_name) {
  class_surface <- surface_results[surface_results$class == class_name, ]
  class_reference <- reference_results[reference_results$class == class_name, ]

  z_matrix <- matrix(
    class_surface$estimate,
    nrow = length(predator_levels),
    ncol = length(predator_levels)
  )
  cell_values <- (
    z_matrix[-nrow(z_matrix), -ncol(z_matrix)] +
      z_matrix[-1, -ncol(z_matrix)] +
      z_matrix[-nrow(z_matrix), -1] +
      z_matrix[-1, -1]
  ) / 4

  graphics::par(mar = c(3.2, 2.4, 3.0, 1.1))
  perspective_matrix <- graphics::persp(
    x = predator_levels,
    y = predator_levels,
    z = z_matrix,
    xlim = c(-1, 1),
    ylim = c(-1, 1),
    zlim = z_limits,
    theta = -38,
    phi = 24,
    expand = 0.78,
    col = as.vector(abundance_colors(cell_values)),
    border = grDevices::adjustcolor("white", alpha.f = 0.28),
    shade = 0.12,
    ticktype = "detailed",
    nticks = 3,
    xlab = "Jackal activity (SD)",
    ylab = "Wolf activity (SD)",
    zlab = "Relative abundance",
    main = class_name,
    cex.axis = 0.72,
    cex.lab = 0.82,
    cex.main = 1.08,
    font.main = 2
  )

  projected <- grDevices::trans3d(
    x = class_reference$jackal_z,
    y = class_reference$wolf_z,
    z = class_reference$estimate,
    pmat = perspective_matrix
  )
  graphics::points(
    projected$x,
    projected$y,
    pch = 21,
    cex = 0.48,
    lwd = 0.55,
    col = "#303030",
    bg = "white"
  )
}

draw_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.16))
  for (class_name in unname(class_labels)) {
    draw_surface_panel(class_name)
  }
  draw_color_key()
}

grDevices::png(
  paste0(output_prefix, ".png"),
  width = 12.4,
  height = 4.9,
  units = "in",
  res = 400,
  type = "cairo-png",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    paste0(output_prefix, ".pdf"),
    width = 12.4,
    height = 4.9,
    family = "sans",
    bg = "white"
  )
} else {
  grDevices::pdf(
    paste0(output_prefix, ".pdf"),
    width = 12.4,
    height = 4.9,
    family = "sans",
    useDingbats = FALSE,
    bg = "white"
  )
}
draw_figure()
grDevices::dev.off()

message("Verified all 27 heatmap reference predictions against: ", reference_grid_file)
message("Wrote: ", paste0(output_prefix, ".png"))
message("Wrote: ", paste0(output_prefix, ".pdf"))
message("Wrote: ", paste0(output_prefix, "_surface.csv"))
