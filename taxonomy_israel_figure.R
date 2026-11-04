# Publication figure: recruitment and current composition of academic
# taxonomists in Israel. All values below come from the supplied summary data.

library(ggplot2)
library(patchwork)

# -----------------------------------------------------------------------------
# Data and derived quantities
# -----------------------------------------------------------------------------

recruitment <- data.frame(
  period = factor(
    c("2004\u20132021", "2021\u20132026"),
    levels = c("2004\u20132021", "2021\u20132026")
  ),
  hires = c(2, 9),
  years = c(17, 5),
  group = c("TAU, before 2021", "The Steinhardt Museum, TAU"),
  stringsAsFactors = FALSE
)
recruitment$rate <- recruitment$hires / recruitment$years
recruitment$bar_label <- c(
  sprintf("%d hires / %d years\n%.2f per year",
          recruitment$hires[1], recruitment$years[1], recruitment$rate[1]),
  sprintf("%d hires / ~%d years\n%.1f per year",
          recruitment$hires[2], recruitment$years[2], recruitment$rate[2])
)

composition <- data.frame(
  category = factor(
    c("TAU, recruited before 2021", "TAU, recruited since 2021",
      "Hebrew University"),
    levels = rev(c("TAU, recruited before 2021", "TAU, recruited since 2021",
                   "Hebrew University"))
  ),
  count = c(3, 9, 1),
  group = c("TAU, before 2021", "The Steinhardt Museum, TAU", "Hebrew University"),
  stringsAsFactors = FALSE
)

# Exact recruitment years and counts supplied for currently active TAU senior
# faculty members. The Hebrew University member is omitted from this table
# because no recruitment year was supplied.
hire_year_counts <- data.frame(
  hire_year = c(2003, 2011, 2019, 2022, 2023, 2024),
  count = c(1, 1, 1, 2, 4, 3),
  group = c(
    rep("TAU, before 2021", 3),
    rep("The Steinhardt Museum, TAU", 3)
  ),
  stringsAsFactors = FALSE
)

# Generate one row per person so the unit chart updates when counts change.
dot_data <- do.call(
  rbind,
  lapply(seq_len(nrow(composition)), function(i) {
    data.frame(
      category = composition$category[i],
      dot = seq_len(composition$count[i]),
      group = composition$group[i]
    )
  })
)
dot_data$category <- factor(dot_data$category, levels = levels(composition$category))

total_taxonomists <- sum(composition$count)
max_dots <- max(composition$count)

# Muted, colour-blind-friendly palette used consistently across panels.
group_colours <- c(
  "TAU, before 2021" = "#4477AA",
  "The Steinhardt Museum, TAU" = "#CC6677",
  "Hebrew University" = "#777777"
)

# -----------------------------------------------------------------------------
# Shared theme
# -----------------------------------------------------------------------------

theme_taxonomy <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.12),
                                margin = margin(b = 8)),
      plot.subtitle = element_text(colour = "grey30", margin = margin(b = 10)),
      plot.caption = element_text(size = rel(0.78), colour = "grey35",
                                  hjust = 0, margin = margin(t = 10)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      plot.margin = margin(12, 16, 12, 12)
    )
}

# -----------------------------------------------------------------------------
# Panel A: period-level recruitment rates at Tel Aviv University
# -----------------------------------------------------------------------------

panel_a <- ggplot(recruitment, aes(x = period, y = rate, fill = group)) +
  geom_col(width = 0.48, show.legend = FALSE) +
  geom_text(
    aes(label = bar_label),
    vjust = -0.35, lineheight = 0.95, size = 3.5, fontface = "bold"
  ) +
  annotate(
    "segment", x = 1.5, xend = 1.5, y = 0.12, yend = 0.84,
    linewidth = 0.45, linetype = "dashed", colour = "grey50"
  ) +
  annotate(
    "text", x = 1.5, y = 0.90,
    label = "2021: Steinhardt Museum recognized\nas an independent academic unit",
    size = 3.05, lineheight = 0.95, colour = "grey25"
  ) +
  annotate(
    "text", x = 1.03, y = 0.51,
    label = paste0("3 hires if the taxonomist recruited just before\n",
                   "the VATAT support model is included"),
    size = 2.75, lineheight = 0.95, colour = "grey40", hjust = 0.5
  ) +
  annotate(
    "text", x = 1.58, y = 2.18,
    label = ">10-fold increase in\nTAU recruitment rate",
    size = 3.5, fontface = "bold", lineheight = 0.95, colour = "grey15"
  ) +
  scale_fill_manual(values = group_colours) +
  scale_y_continuous(
    name = "Active taxonomists recruited per year",
    limits = c(0, 2.38), breaks = seq(0, 2, by = 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "A. Recruitment of academic taxonomists at Tel Aviv University",
    x = NULL
  ) +
  theme_taxonomy() +
  theme(axis.line.x = element_line(colour = "grey45", linewidth = 0.4))

# -----------------------------------------------------------------------------
# Panel B: current senior academic composition (one point per taxonomist)
# -----------------------------------------------------------------------------

panel_b <- ggplot(dot_data, aes(x = dot, y = category, colour = group)) +
  geom_point(size = 7.2, show.legend = FALSE) +
  annotate(
    "text", x = max_dots + 0.55, y = 0.62,
    label = sprintf("Total n = %d", total_taxonomists),
    hjust = 1, size = 3.8, fontface = "bold", colour = "grey20"
  ) +
  scale_colour_manual(values = group_colours) +
  scale_x_continuous(
    limits = c(0.55, max_dots + 0.55),
    breaks = seq_len(max_dots),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B. Current active academic taxonomists",
    subtitle = "One dot = one active senior academic taxonomist",
    x = NULL, y = NULL,
    caption = paste0(
      "Note: An additional active taxonomist holds a research position at IOLR\n",
      "and is not included in the academic total."
    )
  ) +
  theme_taxonomy() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.caption = element_text(
      size = 8.6, colour = "grey35", hjust = 0, margin = margin(t = 10)
    )
  )

# -----------------------------------------------------------------------------
# Assemble and export both requested layouts
# -----------------------------------------------------------------------------

horizontal_figure <- panel_a + panel_b +
  plot_layout(widths = c(1.08, 1)) &
  theme(plot.margin = margin(14, 16, 14, 14))

vertical_figure <- panel_a / panel_b +
  plot_layout(heights = c(1.08, 1)) &
  theme(plot.margin = margin(14, 16, 14, 14))

output_dir <- "figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "taxonomy_israel_panels_horizontal.png"),
       horizontal_figure, width = 12, height = 6.5, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(output_dir, "taxonomy_israel_panels_vertical.png"),
       vertical_figure, width = 8, height = 10, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(output_dir, "taxonomy_israel_panel_A.png"),
       panel_a, width = 6.5, height = 6, units = "in",
       dpi = 300, bg = "white")
ggsave(file.path(output_dir, "taxonomy_israel_panel_B.png"),
       panel_b, width = 6.5, height = 6, units = "in",
       dpi = 300, bg = "white")

# Display both alternatives when run interactively or sourced in an IDE.
print(horizontal_figure)
print(vertical_figure)

# -----------------------------------------------------------------------------
# Supplementary timeline: cumulative current workforce by recruitment cohort
# -----------------------------------------------------------------------------

# Each TAU cohort accumulates from its supplied recruitment year to the present.
# The HU member is shown only at the current endpoint because no hire year was
# supplied for that person.
years_shown <- 2003:2026
hu_current <- composition$count[composition$group == "Hebrew University"]
israel_current <- sum(hire_year_counts$count) + hu_current

# Generate cumulative annual counts directly from the hire-year table.
timeline_bars <- rbind(
  do.call(rbind, lapply(unique(hire_year_counts$group), function(group_name) {
    group_hires <- hire_year_counts[hire_year_counts$group == group_name, ]
    data.frame(
      year = years_shown,
      component = group_name,
      count = vapply(
        years_shown,
        function(y) sum(group_hires$count[group_hires$hire_year <= y]),
        numeric(1)
      )
    )
  })),
  data.frame(year = years_shown, component = "Hebrew University",
             count = ifelse(years_shown == 2026, hu_current, 0))
)
timeline_bars$component <- factor(
  timeline_bars$component,
  levels = c("Hebrew University", "The Steinhardt Museum, TAU",
             "TAU, before 2021")
)
# Place each annual bar within its calendar-year interval rather than centering
# it on the year boundary. This keeps the October 2021 marker before 2022 hires.
timeline_bars$bar_x <- timeline_bars$year + 0.41
bar_width <- 0.718
vatat_marker_x <- mean(c(2003, 2004)) + 0.41
museum_marker_x <- (2021 + 0.41) - bar_width / 2 + 0.75 * bar_width

timeline_figure <- ggplot(
  timeline_bars, aes(x = bar_x, y = count, fill = component)
) +
  geom_col(width = bar_width, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(vatat_marker_x, museum_marker_x), colour = "grey35",
             linewidth = 0.5, linetype = "dashed") +
  annotate("curve", x = 2004.55, y = -1.25,
           xend = vatat_marker_x, yend = 0.12,
           curvature = -0.18, arrow = arrow(length = grid::unit(0.12, "inches")),
           linewidth = 0.45, colour = "grey30") +
  annotate("text", x = 2004.35, y = -1.72,
           label = "VATAT support\nbegins",
           hjust = 0, vjust = 0.5, size = 2.85, lineheight = 0.95,
           colour = "grey20") +
  annotate("curve", x = 2019.35, y = -1.65,
           xend = museum_marker_x, yend = 0.12,
           curvature = -0.18, arrow = arrow(length = grid::unit(0.12, "inches")),
           linewidth = 0.45, colour = "grey30") +
  annotate("text", x = 2019.25, y = -1.82,
           label = "The Steinhardt Museum recognized\nas an independent academic unit",
           hjust = 1, vjust = 0.5, size = 2.85, lineheight = 0.95,
           colour = "grey20") +
  scale_fill_manual(
    values = group_colours, name = NULL,
    breaks = c("TAU, before 2021", "The Steinhardt Museum, TAU",
               "Hebrew University"),
    labels = c("TAU", "The Steinhardt Museum", "The Hebrew University")
  ) +
  scale_x_continuous(
    breaks = c(2006, 2011, 2016, 2021, 2026) + 0.41,
    labels = c(2006, 2011, 2016, 2021, 2026),
    limits = c(2003, 2027),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = c(0, 3, 6, 9, 12, 13), limits = c(-2.45, 14.5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Current senior faculty members in taxonomy\nby recruitment period",
    x = NULL,
    y = "Number of active\nsenior faculty members"
  ) +
  theme_taxonomy(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.line.x = element_line(colour = "grey45", linewidth = 0.4),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 9.5)
  )

ggsave(file.path(output_dir, "taxonomy_israel_cumulative_timeline.png"),
       timeline_figure, width = 5.2, height = 5.2, units = "in",
       dpi = 300, bg = "white")
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
ggsave(file.path(output_dir, "taxonomy_israel_cumulative_timeline.pdf"),
       timeline_figure, width = 5.2, height = 5.2, units = "in",
       device = pdf_device, bg = "white")

print(timeline_figure)
