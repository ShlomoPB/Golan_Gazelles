#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(spAbundance))

parse_args <- function(args) {
  cfg <- list(
    period = "Overall",
    mode = "full",
    data = "gazelle_multispecies_NMix_data.rds",
    obs_covs = "ObsCovs_050524.csv",
    events = "Gazella gazella_independent_events_30min.csv",
    refit = "false",
    seed = "20260804"
  )
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) == 2) cfg[[pieces[1]]] <- pieces[2]
  }
  cfg
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
valid_periods <- c("Overall", "Early", "Late")
if (!cfg$period %in% valid_periods) {
  stop("--period must be Overall, Early, or Late")
}
if (!cfg$mode %in% c("smoke", "full")) {
  stop("--mode must be smoke or full")
}

settings <- if (cfg$mode == "smoke") {
  list(n.batch = 20L, batch.length = 10L, n.burn = 50L, n.thin = 2L, n.chains = 2L)
} else {
  list(n.batch = 2000L, batch.length = 25L, n.burn = 5000L, n.thin = 5L, n.chains = 3L)
}

period_slug <- c(Overall = "overall", Early = "early_hotter", Late = "late_mild_cold")[[cfg$period]]
fit_tag <- if (cfg$mode == "smoke") {
  "smoke"
} else {
  paste0(
    settings$n.batch, "b_", settings$batch.length, "bl_",
    settings$n.thin, "th_", settings$n.chains, "ch"
  )
}
output_file <- paste0(
  "gazelle_lfMsNMix_2f_NB_wolf_jackal_interaction_",
  period_slug, "_", fit_tag, "_fit.rds"
)

if (file.exists(output_file) && tolower(cfg$refit) != "true") {
  cat("Fit already exists; leaving it unchanged:", output_file, "\n")
  quit(save = "no", status = 0)
}

make_det_model_frame <- function(det_covs) {
  out <- as.data.frame(
    lapply(det_covs, function(x) as.vector(x)),
    stringsAsFactors = FALSE
  )
  for (nm in names(out)) {
    if (is.character(out[[nm]])) out[[nm]] <- factor(out[[nm]])
  }
  out
}

build_season_matrix <- function(data_list, obs_cov_file, event_file) {
  events <- read.csv(event_file, check.names = FALSE, stringsAsFactors = FALSE)
  events$ts <- as.POSIXct(
    events$event_time,
    format = "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  )
  split_date <- as.Date(stats::median(events$ts, na.rm = TRUE), tz = "UTC")

  obs <- read.csv(
    obs_cov_file,
    na.strings = c("NA", ""),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  obs$camera_id <- as.integer(obs[["Camera Number"]])
  obs$round <- as.integer(obs[["Round Number"]])
  obs$start_date <- as.Date(obs[["Start date"]], format = "%d/%m/%Y")
  obs$end_date <- as.Date(obs[["End Date"]], format = "%d/%m/%Y")
  obs$mid_date <- obs$start_date + floor(as.numeric(obs$end_date - obs$start_date) / 2)

  round_mid <- stats::aggregate(
    as.numeric(obs$mid_date),
    list(round = obs$round),
    function(x) stats::median(x, na.rm = TRUE)
  )
  names(round_mid)[2] <- "round_mid"
  obs <- merge(obs, round_mid, by = "round", all.x = TRUE, sort = FALSE)
  missing_mid <- is.na(obs$mid_date)
  obs$mid_date[missing_mid] <- as.Date(obs$round_mid[missing_mid], origin = "1970-01-01")
  obs$period <- ifelse(obs$mid_date <= split_date, "Early", "Late")

  site_ids <- as.integer(dimnames(data_list$y)[[2]])
  round_ids <- as.integer(sub("^X", "", dimnames(data_list$y)[[3]]))
  expected <- expand.grid(
    camera_id = site_ids,
    round = round_ids,
    KEEP.OUT.ATTRS = FALSE
  )
  expected$key <- paste(expected$camera_id, expected$round, sep = ":")
  obs$key <- paste(obs$camera_id, obs$round, sep = ":")
  if (anyDuplicated(obs$key)) stop("Observation covariates contain duplicate camera-round records.")
  match_index <- match(expected$key, obs$key)
  if (anyNA(match_index)) stop("Could not assign a season to every camera-round.")

  season_mat <- matrix(
    obs$period[match_index],
    nrow = length(site_ids),
    ncol = length(round_ids),
    byrow = FALSE,
    dimnames = dimnames(data_list$y)[2:3]
  )
  attr(season_mat, "split_date") <- split_date
  season_mat
}

mask_period <- function(data_list, season_mat, period) {
  keep_mask <- season_mat == period
  out <- data_list
  for (s in seq_len(dim(out$y)[1])) {
    y_s <- out$y[s, , ]
    y_s[!keep_mask] <- NA
    out$y[s, , ] <- y_s
  }
  for (nm in names(out$det.covs)) {
    cov_i <- out$det.covs[[nm]]
    cov_i[!keep_mask] <- NA
    out$det.covs[[nm]] <- cov_i
  }
  keep_sites <- apply(!is.na(out$y), 2, any)
  out$y <- out$y[, keep_sites, , drop = FALSE]
  out$abund.covs <- out$abund.covs[keep_sites, , drop = FALSE]
  out$coords <- out$coords[keep_sites, , drop = FALSE]
  for (nm in names(out$det.covs)) {
    out$det.covs[[nm]] <- out$det.covs[[nm]][keep_sites, , drop = FALSE]
  }
  out
}

data_list <- readRDS(cfg$data)
season_mat <- build_season_matrix(data_list, cfg$obs_covs, cfg$events)
if (cfg$period != "Overall") {
  data_list <- mask_period(data_list, season_mat, cfg$period)
}

abund_formula <- ~ Null.Wolf.Abundance.z * null.jackal.abundance.z
det_formula <- ~ path_type + camera_model + vegetation_structure + area_behind_camera
X_abund <- model.matrix(abund_formula, data = data_list$abund.covs)
X_det <- model.matrix(det_formula, data = make_det_model_frame(data_list$det.covs))
n_species <- dim(data_list$y)[1]

priors <- list(
  beta.comm.normal = list(mean = rep(0, ncol(X_abund)), var = rep(100, ncol(X_abund))),
  alpha.comm.normal = list(mean = rep(0, ncol(X_det)), var = rep(2.72, ncol(X_det))),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
  kappa.unif = list(a = rep(0, n_species), b = rep(10, n_species))
)

N_init <- apply(data_list$y, c(1, 2), max, na.rm = TRUE)
N_init[!is.finite(N_init)] <- 0
inits <- list(
  beta.comm = 0,
  alpha.comm = 0,
  beta = 0,
  alpha = 0,
  tau.sq.beta = 0.5,
  tau.sq.alpha = 0.5,
  N = N_init,
  kappa = rep(1, n_species)
)
tuning <- list(beta = 0.5, alpha = 0.5, lambda = 0.5, w = 0.5, kappa = 0.5)

cat("Fitting", cfg$period, "joint wolf-jackal model\n")
cat("  sites:", dim(data_list$y)[2], "\n")
cat("  observed camera-round cells:", sum(!is.na(data_list$y[1, , ])), "\n")
cat("  season split date:", as.character(attr(season_mat, "split_date")), "\n")
cat("  mode:", cfg$mode, "\n")

set.seed(as.integer(cfg$seed) + match(cfg$period, valid_periods))
fit <- lfMsNMix(
  abund.formula = abund_formula,
  det.formula = det_formula,
  data = data_list,
  inits = inits,
  priors = priors,
  tuning = tuning,
  n.factors = 2,
  n.batch = settings$n.batch,
  batch.length = settings$batch.length,
  n.burn = settings$n.burn,
  n.thin = settings$n.thin,
  n.chains = settings$n.chains,
  family = "NB",
  n.omp.threads = 1,
  verbose = cfg$mode == "full",
  n.report = max(1L, floor(settings$n.batch / 10))
)

attr(fit, "analysis_period") <- cfg$period
attr(fit, "season_split_date") <- attr(season_mat, "split_date")
attr(fit, "fit_settings") <- settings
saveRDS(fit, output_file)
cat("Saved:", output_file, "\n")
