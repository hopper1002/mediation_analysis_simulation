# Independent implementation of ASHMED, supplied manuscript, Section 2.2.
# Fixed zero-centred same-scale priors; no access to simulation truth.

ashmed_options <- function(options = list()) {
  defaults <- list(J = 12L, r_min = 0.25, r_max = 8,
                   correlations = c(-0.5, 0, 0.5),
                   initial_pi = c(H00 = 0.85, H10 = 0.1, H01 = 0.1, H11 = 0.05),
                   relative_tolerance = 1e-8, kkt_tolerance = 1e-4,
                   max_iter = 5000L, accelerate = TRUE,
                   fit_subset = NULL, posterior_chunk_size = 5000L)
  unknown <- setdiff(names(options), names(defaults))
  if (length(unknown)) stop("Unknown ASHMED options: ", paste(unknown, collapse = ", "))
  opt <- utils::modifyList(defaults, options)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  integer_scalar <- function(x) scalar(x) && x == floor(x)
  stopifnot(integer_scalar(opt$J), opt$J >= 2L,
            scalar(opt$r_min), scalar(opt$r_max), opt$r_min > 0, opt$r_max > opt$r_min,
            is.numeric(opt$correlations), length(opt$correlations) > 0L,
            all(is.finite(opt$correlations)), all(abs(opt$correlations) < 1),
            !anyDuplicated(opt$correlations),
            scalar(opt$relative_tolerance), opt$relative_tolerance > 0,
            scalar(opt$kkt_tolerance), opt$kkt_tolerance > 0,
            integer_scalar(opt$max_iter), opt$max_iter >= 1L,
            is.logical(opt$accelerate), length(opt$accelerate) == 1L, !is.na(opt$accelerate),
            integer_scalar(opt$posterior_chunk_size), opt$posterior_chunk_size >= 1L)
  states <- c("H00", "H10", "H01", "H11")
  stopifnot(is.numeric(opt$initial_pi), length(opt$initial_pi) == 4L,
            identical(names(opt$initial_pi), states), all(is.finite(opt$initial_pi)),
            all(opt$initial_pi > 0))
  # The manuscript's illustrative initial values sum to 1.10, not 1.
  # This normalization changes initialization only, not the fitted model.
  opt$initial_pi <- opt$initial_pi / sum(opt$initial_pi)
  if (!is.null(opt$fit_subset))
    stopifnot(integer_scalar(opt$fit_subset), opt$fit_subset >= 10L)
  opt
}

ashmed_validate_estimates <- function(e) {
  required <- c("alpha_hat", "beta_hat", "var_alpha", "var_beta")
  if (!is.data.frame(e) || !all(required %in% names(e)) || !nrow(e))
    stop("ASHMED requires a nonempty coefficient/variance data frame.")
  if (any(!is.finite(as.matrix(e[, required]))) || any(e$var_alpha <= 0 | e$var_beta <= 0))
    stop("ASHMED requires finite estimates and positive sampling variances.")
  covariance <- if ("cov_alpha_beta" %in% names(e)) e$cov_alpha_beta else rep(0, nrow(e))
  if (!is.numeric(covariance) || any(!is.finite(covariance)) ||
      any(abs(covariance) >= sqrt(e$var_alpha * e$var_beta)))
    stop("Sampling covariance matrices must be positive definite.")
  covariance
}

ashmed_grid <- function(e, opt) {
  r <- exp(seq(log(opt$r_min), log(opt$r_max), length.out = opt$J))
  a <- (r * stats::median(sqrt(e$var_alpha)))^2
  b <- (r * stats::median(sqrt(e$var_beta)))^2
  component <- function(state, j, rho, u, v, c) {
    data.frame(state = state, scale_index = j, rho = rho,
               prior_var_alpha = u, prior_var_beta = v, prior_cov = c)
  }
  grid <- rbind(component("H00", 0L, 0, 0, 0, 0),
                component("H10", seq_along(r), 0, a, 0, 0),
                component("H01", seq_along(r), 0, 0, b, 0))
  for (rho in opt$correlations)
    grid <- rbind(grid, component("H11", seq_along(r), rho, a, b, rho * sqrt(a * b)))
  grid$component_id <- sprintf("%s_j%02d_rho%+.2f", grid$state, grid$scale_index, grid$rho)
  grid$scale_multiplier <- ifelse(grid$scale_index == 0, 0, r[pmax(1L, grid$scale_index)])
  rownames(grid) <- NULL
  grid
}

ashmed_log_likelihood <- function(e, grid, sampling_cov = ashmed_validate_estimates(e)) {
  p <- nrow(e)
  answer <- matrix(0, p, nrow(grid))
  for (k in seq_len(nrow(grid))) {
    a <- e$var_alpha + grid$prior_var_alpha[k]
    b <- e$var_beta + grid$prior_var_beta[k]
    c <- sampling_cov + grid$prior_cov[k]
    determinant <- a * b - c^2
    if (any(!is.finite(determinant) | determinant <= 0)) stop("Invalid marginal covariance.")
    quadratic <- (b * e$alpha_hat^2 - 2 * c * e$alpha_hat * e$beta_hat + a * e$beta_hat^2) / determinant
    answer[, k] <- -log(2 * pi) - 0.5 * (log(determinant) + quadratic)
  }
  if (any(!is.finite(answer))) stop("Non-finite marginal log likelihood.")
  colnames(answer) <- grid$component_id
  answer
}

ashmed_fit_weights <- function(log_likelihood, grid, opt) {
  # Flatten w_k = pi_state * omega_k. EM then has w_new = mean(r_ik),
  # exactly the product of the manuscript's two M-step updates.
  offset <- apply(log_likelihood, 1L, max)
  F <- exp(log_likelihood - offset)
  p <- nrow(F)
  counts <- table(grid$state)
  w <- unname(opt$initial_pi[grid$state] / counts[grid$state])
  w <- w / sum(w)
  objective <- function(x) {
    den <- as.vector(F %*% x)
    if (any(!is.finite(den) | den <= 0)) return(-Inf)
    sum(log(den)) # Row offsets are constant during weight optimization.
  }
  em <- function(x) {
    den <- as.vector(F %*% x)
    next_w <- x * as.vector(crossprod(F, 1 / den)) / p
    next_w / sum(next_w)
  }
  value <- objective(w)
  trace <- data.frame(iteration = 0L, loglik = value + sum(offset), dual_gap = NA_real_,
                      relative_change = NA_real_, update = "initial")
  converged <- FALSE
  for (iteration in seq_len(opt$max_iter)) {
    w1 <- em(w)
    candidate <- w1
    update <- "EM"
    if (opt$accelerate) {
      w2 <- em(w1)
      candidate <- w2
      update <- "EM2"
      r <- w1 - w
      v <- w2 - w1 - r
      if (sum(v^2) > 0 && sum(r^2) > 0) {
        step <- -max(1, min(1000, sqrt(sum(r^2) / sum(v^2))))
        base_value <- objective(w2)
        for (attempt in seq_len(15L)) {
          extrapolated <- w - 2 * step * r + step^2 * v
          # Projection is a numerical proposal, accepted only if the same
          # unpenalized likelihood improves over the ordinary EM2 update.
          extrapolated <- pmax(0, extrapolated)
          extrapolated <- extrapolated / sum(extrapolated)
          if (is.finite(objective(extrapolated))) {
            proposed <- em(extrapolated)
            if (objective(proposed) >= base_value) {
              candidate <- proposed
              update <- "accelerated_EM"
              break
            }
          }
          step <- (step - 1) / 2
        }
      }
    }
    next_value <- objective(candidate)
    denominator <- as.vector(F %*% candidate)
    gradient <- as.vector(crossprod(F, 1 / denominator)) / p
    dual_gap <- max(0, max(gradient) - 1)
    # Fixed-grid log likelihood is concave in w. A vertex line search can
    # revive a component projected to zero and resolve stalled EM updates.
    if (dual_gap > opt$kkt_tolerance &&
        ((next_value - value) / max(1, abs(value + sum(offset))) < opt$relative_tolerance || iteration %% 20L == 0L)) {
      best <- which.max(gradient)
      direction <- F[, best] - denominator
      line <- stats::optimize(function(t) -sum(log(denominator + t * direction)),
                              interval = c(0, 1 - 1e-12), tol = 1e-10)
      if (-line$objective > next_value) {
        candidate <- (1 - line$minimum) * candidate
        candidate[best] <- candidate[best] + line$minimum
        next_value <- -line$objective
        denominator <- as.vector(F %*% candidate)
        gradient <- as.vector(crossprod(F, 1 / denominator)) / p
        dual_gap <- max(0, max(gradient) - 1)
        update <- paste0(update, "+vertex")
      }
    }
    if (next_value < value - 1e-9 * max(1, abs(value))) stop("ASHMED likelihood decreased.")
    relative <- abs(next_value - value) / max(1, abs(value + sum(offset)))
    w <- candidate / sum(candidate)
    value <- next_value
    trace[iteration + 1L, ] <- list(iteration, value + sum(offset), dual_gap, relative, update)
    if (relative <= opt$relative_tolerance && dual_gap <= opt$kkt_tolerance) {
      converged <- TRUE
      break
    }
  }
  list(weights = w, loglik = value + sum(offset), trace = trace,
       iterations = iteration, converged = converged, dual_gap = dual_gap)
}

ashmed_responsibilities <- function(log_likelihood, weights) {
  offset <- apply(log_likelihood, 1L, max)
  values <- sweep(exp(log_likelihood - offset), 2L, weights, "*")
  denominator <- rowSums(values)
  if (any(!is.finite(denominator) | denominator <= 0)) stop("Invalid posterior normalizer.")
  values / denominator
}

ashmed_posterior_moments <- function(e, grid, responsibilities, sampling_cov) {
  p <- nrow(e)
  mean_a <- mean_b <- second_a <- second_b <- mean_ab <- second_ab <- numeric(p)
  for (k in seq_len(nrow(grid))) {
    u <- grid$prior_var_alpha[k]
    v <- grid$prior_var_beta[k]
    c <- grid$prior_cov[k]
    if (u == 0 && v == 0) next
    a <- e$var_alpha + u
    b <- e$var_beta + v
    total_c <- sampling_cov + c
    determinant <- a * b - total_c^2
    # Conjugate Gaussian formulas use (S + U)^-1, never U^-1;
    # valid for axis components and the point mass as well.
    B11 <- (u * b - c * total_c) / determinant
    B12 <- (-u * total_c + c * a) / determinant
    B21 <- (c * b - v * total_c) / determinant
    B22 <- (-c * total_c + v * a) / determinant
    ma <- B11 * e$alpha_hat + B12 * e$beta_hat
    mb <- B21 * e$alpha_hat + B22 * e$beta_hat
    Vaa <- pmax(0, u - B11 * u - B12 * c)
    Vbb <- pmax(0, v - B21 * c - B22 * v)
    Vab <- ((c - B11 * c - B12 * v) + (c - B21 * u - B22 * c)) / 2
    rk <- responsibilities[, k]
    mean_a <- mean_a + rk * ma
    mean_b <- mean_b + rk * mb
    second_a <- second_a + rk * (Vaa + ma^2)
    second_b <- second_b + rk * (Vbb + mb^2)
    mean_ab <- mean_ab + rk * (Vab + ma * mb)
    second_ab <- second_ab + rk * (ma^2 * mb^2 + ma^2 * Vbb + mb^2 * Vaa +
                                       4 * ma * mb * Vab + Vaa * Vbb + 2 * Vab^2)
  }
  data.frame(alpha_mean = mean_a, alpha_sd = sqrt(pmax(0, second_a - mean_a^2)),
             beta_mean = mean_b, beta_sd = sqrt(pmax(0, second_b - mean_b^2)),
             alpha_beta_cov = mean_ab - mean_a * mean_b,
             mediation_mean = mean_ab,
             mediation_sd = sqrt(pmax(0, second_ab - mean_ab^2)))
}

fit_ashmed <- function(e, options = list()) {
  opt <- ashmed_options(options)
  sampling_cov <- ashmed_validate_estimates(e)
  # Grid uses median SE of all observed tests even when fitting a subset.
  grid <- ashmed_grid(e, opt)
  fit_indices <- seq_len(nrow(e))
  if (!is.null(opt$fit_subset) && opt$fit_subset < nrow(e))
    fit_indices <- sort(sample.int(nrow(e), opt$fit_subset))
  fit <- ashmed_fit_weights(ashmed_log_likelihood(e[fit_indices, , drop = FALSE], grid,
                                                sampling_cov[fit_indices]), grid, opt)
  states <- c("H00", "H10", "H01", "H11")
  pi_state <- setNames(vapply(states, function(s) sum(fit$weights[grid$state == s]), numeric(1)), states)
  omega <- fit$weights
  for (s in states) {
    idx <- grid$state == s
    omega[idx] <- if (pi_state[s] > 0) fit$weights[idx] / pi_state[s] else 1 / sum(idx)
  }
  grid$weight <- fit$weights
  grid$within_state_weight <- omega
  chunks <- split(seq_len(nrow(e)), ceiling(seq_len(nrow(e)) / opt$posterior_chunk_size))
  post <- effects <- list()
  for (j in seq_along(chunks)) {
    idx <- chunks[[j]]
    e_chunk <- e[idx, , drop = FALSE]
    resp <- ashmed_responsibilities(ashmed_log_likelihood(e_chunk, grid, sampling_cov[idx]), fit$weights)
    post[[j]] <- vapply(states, function(s) rowSums(resp[, grid$state == s, drop = FALSE]), numeric(length(idx)))
    # vapply returns a vector when length(idx)==1; restore the p x 4 shape.
    post[[j]] <- matrix(post[[j]], nrow = length(idx), dimnames = list(NULL, states))
    effects[[j]] <- ashmed_posterior_moments(e_chunk, grid, resp, sampling_cov[idx])
  }
  posterior <- do.call(rbind, post)
  moments <- do.call(rbind, effects)
  rownames(moments) <- NULL
  list(lfdr = pmin(1, pmax(0, rowSums(posterior[, c("H00", "H10", "H01"), drop = FALSE]))),
       posterior = posterior, moments = moments, components = grid, pi = pi_state,
       fit = fit, options = opt, fit_indices = fit_indices,
       sampling_covariance = if (any(sampling_cov != 0)) "provided_full" else "working_diagonal")
}

run_ashmed <- function(input, q, options = list()) {
  fit <- fit_ashmed(input$estimates, options)
  if (!fit$fit$converged)
    warning(sprintf("ASHMED optimizer reached max_iter; dual gap %.3g. Inspect fit diagnostics.", fit$fit$dual_gap))
  reject <- lfdr_step_up(fit$lfdr, q)
  pathway_id <- if ("pathway_id" %in% names(input$estimates)) input$estimates$pathway_id else seq_along(reject)
  posterior_table <- data.frame(pathway_id = pathway_id, fit$posterior, lfdr = fit$lfdr,
                                reject = reject, fit$moments, check.names = FALSE)
  list(reject = reject, score = fit$lfdr, score_type = "local_FDR",
       posterior = posterior_table, components = fit$components,
       diagnostics = list(engine = "ashmed_fixed_grid", pi = fit$pi,
                          converged = fit$fit$converged, iterations = fit$fit$iterations,
                          loglik = fit$fit$loglik, dual_gap = fit$fit$dual_gap,
                          trace = fit$fit$trace, options = fit$options,
                          n_fit = length(fit$fit_indices), n_total = length(reject),
                          sampling_covariance = fit$sampling_covariance,
                          likelihood = if (identical(input$meta$scenario, "binary"))
                            "logistic_summary_normal_extension" else "linear_summary_normal",
                          selected_mean_lfdr = if (any(reject)) mean(fit$lfdr[reject]) else NA_real_))
}

register_plugin <- function(registry) {
  register_method(registry, "ASHMED", run_ashmed, version = "1.0.0",
                  description = "Independent ASHMED Section 2.2: fixed zero-centred four-state shrinkage mixture")
}
