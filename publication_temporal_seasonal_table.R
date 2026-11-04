if (!requireNamespace("readr", quietly = TRUE)) {
  stop("Package 'readr' is required.")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

library(readr)
library(dplyr)

fmt_num <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  paste0(fmt_num(mean, digits), " [", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), "]")
}

overall <- read_csv("gazelle_temporal_activity_solar_season_overall_class_comparisons.csv", show_col_types = FALSE)
hl <- read_csv("gazelle_temporal_activity_solar_season_high_low_pressure_comparisons.csv", show_col_types = FALSE)

overall_tbl <- overall %>%
  transmute(
    period,
    comparison,
    n1,
    n2,
    overlap = ifelse(is.na(overlap), NA_character_, fmt_ci(overlap, overlap_lwr, overlap_upr)),
    activity_1 = ifelse(is.na(activity1), NA_character_, fmt_num(activity1)),
    activity_2 = ifelse(is.na(activity2), NA_character_, fmt_num(activity2)),
    p_value = ifelse(is.na(activity_diff_p), NA_character_, fmt_num(activity_diff_p))
  )

hl_tbl <- hl %>%
  filter(!is.na(activity_diff_p), activity_diff_p <= 0.05) %>%
  arrange(period, class, overlap, activity_diff_p) %>%
  group_by(period, class) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  transmute(
    period,
    class,
    covariate = covariate_label,
    low_events = n1,
    high_events = n2,
    overlap = fmt_ci(overlap, overlap_lwr, overlap_upr),
    p_value = fmt_num(activity_diff_p)
  )

write_csv(overall_tbl, "publication_temporal_seasonal_overall_table.csv")
write_csv(hl_tbl, "publication_temporal_seasonal_shift_table.csv")

lines <- c(
  "# Seasonal Temporal Activity Tables",
  "",
  "## Table T1. Overall class comparisons in apparent solar time",
  "",
  "| Period | Comparison | n1 | n2 | Overlap [95% CI] | Activity 1 | Activity 2 | Kernel comparison p |",
  "|---|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(overall_tbl))) {
  lines <- c(
    lines,
    paste0(
      "| `", overall_tbl$period[i], "` | `", overall_tbl$comparison[i], "` | `", overall_tbl$n1[i],
      "` | `", overall_tbl$n2[i], "` | `", overall_tbl$overlap[i], "` | `",
      overall_tbl$activity_1[i], "` | `", overall_tbl$activity_2[i], "` | `", overall_tbl$p_value[i], "` |"
    )
  )
}

lines <- c(
  lines,
  "",
  "## Table T2. Strongest significant seasonal high-vs-low quartile temporal shifts",
  "",
  "| Period | Class | Covariate | Low events | High events | Overlap [95% CI] | Kernel comparison p |",
  "|---|---|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(hl_tbl))) {
  lines <- c(
    lines,
    paste0(
      "| `", hl_tbl$period[i], "` | `", hl_tbl$class[i], "` | `", hl_tbl$covariate[i], "` | `",
      hl_tbl$low_events[i], "` | `", hl_tbl$high_events[i], "` | `", hl_tbl$overlap[i],
      "` | `", hl_tbl$p_value[i], "` |"
    )
  )
}

writeLines(lines, "publication_temporal_seasonal_table.md")
cat("Saved seasonal temporal tables\n")
