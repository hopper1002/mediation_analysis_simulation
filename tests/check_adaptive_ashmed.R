PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source("R/bootstrap.R", encoding = "UTF-8")
registry <- load_method_registry(PROJECT_ROOT)
env <- environment(registry$ASHMED_Adaptive$run)
eq <- function(a, b, tol = 1e-9) stopifnot(isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE)))
e <- data.frame(alpha_hat = c(-.4, .8), beta_hat = c(.9, -.3),
                var_alpha = c(.1, .15), var_beta = c(.2, .12), cov_alpha_beta = c(.02, -.015))
opt <- env$adaptive_options(list(zero_anchor = TRUE, null_pseudocount = 30))
loc <- list(alpha = list(mean = .3, enabled = TRUE), beta = list(mean = -.7, enabled = TRUE))
grid <- env$adaptive_grid(e, opt, loc)
stopifnot(nrow(grid) == 196L, all(grid$mean_beta[grid$state == "H10"] == 0),
          all(grid$mean_alpha[grid$state == "H01"] == 0))
logL <- env$adaptive_log_likelihood(e, grid)
for (k in c(1, 8, 14, 47, 196)) for (i in 1:2) {
  S <- matrix(c(e$var_alpha[i], e$cov_alpha_beta[i], e$cov_alpha_beta[i], e$var_beta[i]), 2)
  U <- matrix(c(grid$prior_var_alpha[k], grid$prior_cov[k], grid$prior_cov[k], grid$prior_var_beta[k]), 2)
  mu <- c(grid$mean_alpha[k], grid$mean_beta[k])
  y <- c(e$alpha_hat[i], e$beta_hat[i])
  eq(logL[i, k], -log(2 * pi) - .5 * log(det(S + U)) - .5 * drop(t(y - mu) %*% solve(S + U, y - mu)))
  resp <- matrix(0, 2, nrow(grid)); resp[, k] <- 1
  moments <- env$adaptive_moments(e, grid, resp, e$cov_alpha_beta)
  mean <- mu + drop(U %*% solve(S + U, y - mu))
  variance <- U - U %*% solve(S + U, U)
  eq(c(moments$alpha_mean[i], moments$beta_mean[i]), mean)
  eq(c(moments$alpha_sd[i]^2, moments$beta_sd[i]^2), diag(variance))
  eq(moments$mediation_mean[i], prod(mean) + variance[1, 2])
}
# With identical component densities, the MAP state prior has known optimum:
# equal mass on the three null states and zero mass on H11.
toy <- matrix(0, 10, 4)
toygrid <- data.frame(state = c("H00", "H10", "H01", "H11"))
fit <- env$adaptive_fit_weights(toy, toygrid, opt)
stopifnot(fit$converged, fit$weights[4] < 1e-7, all(diff(fit$trace$objective) >= -1e-9))
eq(fit$weights[1:3], rep(1/3, 3), tol = 1e-6)
# Exact equivalence to original ASHMED when locations/penalty are disabled
# and the same-scale original grid settings are used.
manifest <- read_manifest("data/data_c054ed713a76/manifest.rds")
idx <- which(manifest$design$scenario == "linear" & manifest$design$mixture == "dense" & manifest$design$tau == 1.5)[1]
data <- read_dataset(manifest, manifest$design$case_id[idx])
estimates <- estimate_coefficients(data)$estimates
original <- env$ASH_BASE$fit_ashmed(estimates)
same <- env$fit_adaptive_ashmed(estimates, list(learn_location = FALSE, zero_anchor = TRUE, grid_mode = "same", J = 12L,
                correlations = c(-.5, 0, .5), null_pseudocount = 0))
eq(original$lfdr, same$lfdr, tol = 1e-5)
out <- registry$ASHMED_Adaptive$run(make_method_input(data, estimates), .05, list(zero_anchor = FALSE))
invisible(validate_method_result(out, ncol(data$M)))
stopifnot(out$diagnostics$converged, max(abs(rowSums(out$posterior[, c("H00", "H10", "H01", "H11")]) - 1)) < 1e-12,
          !any(out$reject) || mean(out$score[out$reject]) <= .05 + 1e-12)
null_e <- data.frame(alpha_hat = rep(0, 20), beta_hat = rep(0, 20),
                     var_alpha = rep(.1, 20), var_beta = rep(.2, 20))
fallback <- env$fit_adaptive_ashmed(null_e)
stopifnot(fallback$fallback, !any(lfdr_step_up(fallback$lfdr, .05)))
cat("Adaptive ASHMED mathematical checks passed.\n")
