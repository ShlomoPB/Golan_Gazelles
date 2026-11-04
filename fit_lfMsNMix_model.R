library(spAbundance)

# Load prepared data
data_list <- readRDS("gazelle_msNMix_data_list.rds")

# Use Null.Wolf.Abundance as the wolf predictor for the lfMsNMix model.
abund_formula <- ~ Null.Wolf.Abundance.z +
  null.jackal.abundance.z +
  null.boar.abundance.z +
  human.abundance.z +
  Nature.Reserves.z +
  Minefields.z +
  Built.and.Roads.z +
  Firing.Zones.z +
  Proximity.to.Water.z

det_formula <- ~ log_days_operational.z + cattle_sightings.z

p_abund <- ncol(model.matrix(abund_formula, data = data_list$abund.covs))
p_det <- length(attr(terms(det_formula), "term.labels")) + 1

priors <- list(
  beta.comm.normal = list(mean = rep(0, p_abund), var = rep(100, p_abund)),
  alpha.comm.normal = list(mean = rep(0, p_det), var = rep(2.72, p_det)),
  tau.sq.beta.ig = list(a = 0.1, b = 0.1),
  tau.sq.alpha.ig = list(a = 0.1, b = 0.1)
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
  N = N_init
)

tuning <- list(beta = 0.5, alpha = 0.5, kappa = 0.5, w = 0.5, lambda = 0.5)

set.seed(20260511)
cat("Fitting lfMsNMix model with Null.Wolf.Abundance.z and 1 latent factor...\n")

out <- lfMsNMix(
  abund.formula = abund_formula,
  det.formula = det_formula,
  data = data_list,
  inits = inits,
  priors = priors,
  tuning = tuning,
  n.factors = 1,
  family = "NB",
  n.batch = 2000,
  batch.length = 25,
  n.burn = 10000,
  n.thin = 10,
  n.chains = 3,
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 100
)

saveRDS(out, "gazelle_lfMsNMix_fit_nb_null_wolf.rds")

sink("gazelle_lfMsNMix_summary.txt")
print(summary(out, level = "community"))
print(summary(out, level = "species"))
sink()

cat("lfMsNMix model fit complete. Output saved to gazelle_lfMsNMix_fit_nb_null_wolf.rds and summary to gazelle_lfMsNMix_summary.txt\n")