library(spAbundance)

data_list <- readRDS("gazelle_msNMix_data_list.rds")

abund_formula <- ~ Wolf.pack.abundance.z +
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
  kappa = 0.5,
  N = N_init
)

tuning <- list(beta = 0.5, alpha = 0.5, kappa = 0.5)

set.seed(20260510)

# Fit NB model
out_nb <- msNMix(
  abund.formula = abund_formula,
  det.formula = det_formula,
  data = data_list,
  inits = inits,
  priors = priors,
  tuning = tuning,
  family = "NB",  # Negative Binomial
  n.batch = 200,
  batch.length = 25,
  n.burn = 1000,
  n.thin = 1,
  n.chains = 1,  # For quick comparison
  n.omp.threads = 1,
  verbose = TRUE,
  n.report = 50
)

saveRDS(out_nb, "gazelle_msNMix_fit_nb.rds")

# Compute WAIC for both models
fit_poisson <- readRDS("gazelle_msNMix_fit.rds")
waic_poisson <- waicAbund(fit_poisson)
waic_nb <- waicAbund(out_nb)

print("WAIC Poisson:")
print(waic_poisson)
print("WAIC NB:")
print(waic_nb)

# Also PPC for NB
ppc_nb <- ppcAbund(out_nb, fit.stat = 'chi-squared', group = 1)
print("PPC NB:")
print(summary(ppc_nb))