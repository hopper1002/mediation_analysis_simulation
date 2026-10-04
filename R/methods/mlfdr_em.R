# Four-component Gaussian mixture in paper equations (14), (18), (30).
# State order: H00, H10, H01, H11. Only observations enter this fitter.
# Vectorized EM, log densities, multiple deterministic starts, and zero-inclusive
# prior variance bounds implement the paper model without the package's bounds.
mlfdr_log_densities <- function(a, b, va, vb, mu, theta, kappa, psi, pi) {
  la0 <- stats::dnorm(a, 0, sqrt(va), log = TRUE)
  la1 <- stats::dnorm(a, mu, sqrt(va + kappa), log = TRUE)
  lb0 <- stats::dnorm(b, 0, sqrt(vb), log = TRUE)
  lb1 <- stats::dnorm(b, theta, sqrt(vb + psi), log = TRUE)
  sweep(cbind(la0 + lb0, la1 + lb0, la0 + lb1, la1 + lb1), 2, log(pi), "+")
}

mlfdr_posterior <- function(log_density) {
  shift <- apply(log_density, 1, max)
  x <- exp(log_density - shift)
  denom <- rowSums(x)
  list(z = x / denom, loglik = sum(shift + log(denom)))
}

mlfdr_variance_update <- function(x, variance, weights, mean, upper, tolerance) {
  objective <- function(v) sum(weights * stats::dnorm(x, mean, sqrt(variance + v), log = TRUE))
  interior <- stats::optimize(objective, c(0, upper), maximum = TRUE, tol = tolerance)$maximum
  candidates <- c(0, interior, upper)
  candidates[which.max(vapply(candidates, objective, numeric(1)))]
}

fit_mlfdr_paper <- function(alpha, beta, var_alpha, var_beta, n, options = list()) {
  defaults <- list(eps = 1e-4, max_iter = 2000L, n_starts = 3L,
                   variance_tolerance = 1e-6, scale = "sqrt_n")
  opts <- utils::modifyList(defaults, options)
  if (!all(c(length(alpha), length(beta), length(var_alpha), length(var_beta)) == length(alpha)) ||
      !all(is.finite(c(alpha, beta, var_alpha, var_beta))) ||
      any(var_alpha <= 0 | var_beta <= 0)) stop("Invalid input for paper MLFDR EM.")
  stopifnot(opts$eps > 0, opts$max_iter >= 1, opts$n_starts %in% 1:3,
            opts$variance_tolerance > 0, opts$scale %in% c("sqrt_n", "raw"))
  s <- if (opts$scale == "sqrt_n") sqrt(n) else 1
  # Orientation affects initialization only; fitted means remain unrestricted.
  signs <- c(if (mean(alpha) < 0) -1 else 1, if (mean(beta) < 0) -1 else 1)
  a <- alpha * s * signs[1]
  b <- beta * s * signs[2]
  va <- var_alpha * s^2
  vb <- var_beta * s^2
  ua <- max(stats::var(a), max(va), 1 / n * s^2) * 10
  ub <- max(stats::var(b), max(vb), 1 / n * s^2) * 10
  initial_weights <- list(c(.7, .1, .1, .1), c(.4, .2, .2, .2), c(.9, .035, .035, .03))
  starts <- vector("list", opts$n_starts)
  for (start in seq_len(opts$n_starts)) {
    pi <- initial_weights[[start]]
    mu <- unname(stats::quantile(a, .99))
    theta <- unname(stats::quantile(b, .99))
    kappa <- psi <- s^2 / n
    state <- mlfdr_posterior(mlfdr_log_densities(a, b, va, vb, mu, theta, kappa, psi, pi))
    history <- state$loglik
    converged <- FALSE
    for (iteration in seq_len(opts$max_iter)) {
      z <- state$z
      pi_new <- pmax(colMeans(z), 1e-10)
      pi_new <- pi_new / sum(pi_new)
      wa <- z[, 2] + z[, 4]
      wb <- z[, 3] + z[, 4]
      mu_new <- sum(a * wa / (va + kappa)) / sum(wa / (va + kappa))
      theta_new <- sum(b * wb / (vb + psi)) / sum(wb / (vb + psi))
      kappa_new <- mlfdr_variance_update(a, va, wa, mu_new, ua, opts$variance_tolerance)
      psi_new <- mlfdr_variance_update(b, vb, wb, theta_new, ub, opts$variance_tolerance)
      next_state <- mlfdr_posterior(mlfdr_log_densities(a, b, va, vb, mu_new, theta_new, kappa_new, psi_new, pi_new))
      increment <- next_state$loglik - state$loglik
      if (increment < -1e-6) stop("Paper MLFDR EM log likelihood decreased.")
      pi <- pi_new
      mu <- mu_new
      theta <- theta_new
      kappa <- kappa_new
      psi <- psi_new
      state <- next_state
      history <- c(history, state$loglik)
      if (increment <= opts$eps) {
        converged <- TRUE
        break
      }
    }
    starts[[start]] <- list(pi = pi, mu = mu * signs[1] / s, theta = theta * signs[2] / s,
                            kappa = kappa / s^2, psi = psi / s^2,
                            loglik = state$loglik, z = state$z, history = history,
                            iterations = iteration, converged = converged,
                            upper = c(kappa = ua, psi = ub) / s^2)
  }
  # Select using observed likelihood, never simulation labels or FDR/power.
  loglik <- vapply(starts, function(x) x$loglik, numeric(1))
  selected <- which.max(loglik)
  fit <- starts[[selected]]
  if (!fit$converged) warning("Paper MLFDR EM reached max_iter; inspect diagnostics.")
  fit$lfdr <- rowSums(fit$z[, 1:3, drop = FALSE])
  fit$start_loglik <- loglik
  fit$selected_start <- selected
  fit$scale <- opts$scale
  fit$start_converged <- vapply(starts, function(x) x$converged, logical(1))
  fit
}
