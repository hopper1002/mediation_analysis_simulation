# Mathematical and interface checks; these do not generate experiment datasets.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source(file.path(PROJECT_ROOT, "R/bootstrap.R"), encoding = "UTF-8")
registry <- load_method_registry(PROJECT_ROOT)
env <- environment(registry$ASHMED$run)
eq <- function(x, y, tol = 1e-9) stopifnot(isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE)))

e <- data.frame(alpha_hat = c(-1.2, 0, 0.7), beta_hat = c(0.6, 0, -0.8),
                var_alpha = c(0.09, 0.04, 0.16), var_beta = c(0.04, 0.09, 0.25),
                cov_alpha_beta = c(0.01, -0.01, 0.02))
opt <- env$ashmed_options()
grid <- env$ashmed_grid(e, opt)
stopifnot(nrow(grid) == 61L, identical(as.integer(table(grid$state)), c(1L, 12L, 12L, 36L)))
eq(sum(opt$initial_pi), 1)
eq(opt$initial_pi, c(.85, .1, .1, .05) / 1.1)
logL <- env$ashmed_log_likelihood(e, grid)
# Independent generic matrix-density calculation, including nonzero sampling covariance.
for (i in seq_len(nrow(e))) for (k in seq_len(nrow(grid))) {
  S <- matrix(c(e$var_alpha[i], e$cov_alpha_beta[i], e$cov_alpha_beta[i], e$var_beta[i]), 2)
  U <- matrix(c(grid$prior_var_alpha[k], grid$prior_cov[k], grid$prior_cov[k], grid$prior_var_beta[k]), 2)
  theta <- c(e$alpha_hat[i], e$beta_hat[i])
  independent <- -log(2 * pi) - .5 * log(det(S + U)) - .5 * drop(t(theta) %*% solve(S + U, theta))
  eq(logL[i, k], independent)
}

# Posterior moment formulas checked independently for each individual component.
for (k in c(1L, 7L, 16L, 30L, 55L)) {
  responsibility <- matrix(0, nrow(e), nrow(grid))
  responsibility[, k] <- 1
  result <- env$ashmed_posterior_moments(e, grid, responsibility, e$cov_alpha_beta)
  for (i in seq_len(nrow(e))) {
    S <- matrix(c(e$var_alpha[i], e$cov_alpha_beta[i], e$cov_alpha_beta[i], e$var_beta[i]), 2)
    U <- matrix(c(grid$prior_var_alpha[k], grid$prior_cov[k], grid$prior_cov[k], grid$prior_var_beta[k]), 2)
    theta <- c(e$alpha_hat[i], e$beta_hat[i])
    mean <- drop(U %*% solve(S + U, theta))
    variance <- U - U %*% solve(S + U, U)
    eq(c(result$alpha_mean[i], result$beta_mean[i]), mean)
    eq(c(result$alpha_sd[i]^2, result$beta_sd[i]^2), diag(variance))
    eq(result$mediation_mean[i], prod(mean) + variance[1, 2])
    eq(result$mediation_sd[i]^2,
       mean[1]^2 * variance[2, 2] + mean[2]^2 * variance[1, 1] +
         2 * prod(mean) * variance[1, 2] + prod(diag(variance)) + variance[1, 2]^2)
  }
}

# A mixture with an analytically known optimum w=(1/2,1/2).
toy_grid <- data.frame(state = c("H00", "H11"))
toy <- log(matrix(c(.9, .1, .1, .9), nrow = 2))
toy_opt <- env$ashmed_options(list(kkt_tolerance = 1e-8, max_iter = 2000L))
toy_fit <- env$ashmed_fit_weights(toy, toy_grid, toy_opt)
stopifnot(toy_fit$converged, toy_fit$dual_gap < 1e-8)
eq(toy_fit$weights, c(.5, .5), tol = 1e-7)
stopifnot(all(diff(toy_fit$trace$loglik) >= -1e-9))
# Boundary optimum: all mass on one component, without division by a zero state mass.
boundary <- env$ashmed_fit_weights(log(matrix(c(.9, .8, .1, .2), 2)), toy_grid, toy_opt)
stopifnot(boundary$converged, boundary$weights[1] > 1 - 1e-7)

# Real saved coefficients: acceleration and ordinary EM target the same objective.
manifest <- read_manifest(file.path(PROJECT_ROOT, "data/data_c054ed713a76/manifest.rds"))
idx <- which(manifest$design$scenario == "linear" & manifest$design$mixture == "dense" & manifest$design$tau == 1.5)[1]
data <- read_dataset(manifest, manifest$design$case_id[idx])
estimates <- estimate_coefficients(data)$estimates
timing <- system.time(fit <- env$fit_ashmed(estimates))
stopifnot(fit$fit$converged, fit$fit$dual_gap <= fit$options$kkt_tolerance,
          max(abs(rowSums(fit$posterior) - 1)) < 1e-12,
          all(fit$lfdr >= 0 & fit$lfdr <= 1), all(is.finite(as.matrix(fit$moments))),
          all(fit$moments$mediation_sd >= 0), all(diff(fit$fit$trace$loglik) >= -1e-8))
plain <- env$fit_ashmed(estimates, list(accelerate = FALSE, max_iter = 6000L))
stopifnot(plain$fit$converged)
eq(fit$fit$loglik, plain$fit$loglik, tol = 1e-5)
# Row permutation and posterior chunking must preserve path-wise inference.
permutation <- rev(seq_len(nrow(estimates)))
permuted <- env$fit_ashmed(estimates[permutation, ], list(posterior_chunk_size = 333L))
eq(fit$lfdr, permuted$lfdr[permutation], tol = 1e-7)
# A one-row posterior chunk exercises the R simplification edge case.
one <- env$fit_ashmed(estimates[1:10, ], list(posterior_chunk_size = 3L))
stopifnot(identical(dim(one$posterior), c(10L, 4L)))
bad <- e
bad$cov_alpha_beta[1] <- 1
stopifnot(inherits(try(env$fit_ashmed(bad), silent = TRUE), "try-error"))
input <- make_method_input(data, estimates)
stopifnot(!("truth" %in% names(input)))
out <- validate_method_result(registry$ASHMED$run(input, .05), nrow(estimates))
stopifnot(!any(out$reject) || mean(out$score[out$reject]) <= .05 + 1e-12)
cat(sprintf("ASHMED checks passed; saved-data fit: %d iterations, %.3f sec, dual gap %.3g.\n",
            fit$fit$iterations, timing[["elapsed"]], fit$fit$dual_gap))
