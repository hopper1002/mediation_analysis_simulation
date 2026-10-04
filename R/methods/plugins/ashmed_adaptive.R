# ASHMED extension: learned locations, independent scale grids, null-state MAP.
# Original ASHMED is loaded privately and remains unchanged.
ASH_BASE_PATH <- file.path(PROJECT_ROOT, "R/methods/plugins/ashmed.R")
ASH_BASE <- new.env(parent = environment())
sys.source(ASH_BASE_PATH, envir = ASH_BASE)

adaptive_options <- function(options = list()) {
  defaults <- list(J = 6L, r_min = .25, r_max = 8, correlations = 0,
                   grid_mode = "rectangular", learn_location = TRUE,
                   zero_anchor = FALSE, include_atom = TRUE,
                   minimum_location_se = .5, null_pseudocount = 30,
                   relative_tolerance = 1e-8, kkt_tolerance = 1e-4,
                   max_iter = 5000L, accelerate = TRUE,
                   location_max_iter = 1500L, location_tolerance = 1e-6)
  if (length(setdiff(names(options), names(defaults)))) stop("Unknown ASHMED_Adaptive options.")
  opt <- utils::modifyList(defaults, options)
  numeric_scalar <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  stopifnot(opt$grid_mode %in% c("same", "rectangular"),
            all(vapply(opt[c("learn_location", "zero_anchor", "include_atom")],
                       function(x) is.logical(x) && length(x) == 1L && !is.na(x), logical(1))),
            opt$learn_location || opt$zero_anchor,
            numeric_scalar(opt$minimum_location_se), opt$minimum_location_se > 0,
            numeric_scalar(opt$null_pseudocount), opt$null_pseudocount >= 0,
            numeric_scalar(opt$location_max_iter), opt$location_max_iter >= 1,
            opt$location_max_iter == floor(opt$location_max_iter),
            numeric_scalar(opt$location_tolerance), opt$location_tolerance > 0)
  base_keys <- c("J", "r_min", "r_max", "correlations", "relative_tolerance", "kkt_tolerance", "max_iter", "accelerate")
  opt$base <- ASH_BASE$ashmed_options(opt[base_keys])
  opt
}

adaptive_location <- function(x, variance, opt) {
  # Known zero spike plus a free Gaussian slab, used ONLY to learn a basis
  # location. No mediation decisions or latent truth enter this operation.
  scale <- median(sqrt(variance))
  x <- x / scale
  variance <- variance / scale^2
  upper <- 10 * max(stats::var(x), max(variance), 1)
  l0 <- stats::dnorm(x, 0, sqrt(variance), log = TRUE)
  posterior <- function(mu, tau2, pi1) {
    l1 <- stats::dnorm(x, mu, sqrt(variance + tau2), log = TRUE)
    a <- l0 + log1p(-pi1)
    b <- l1 + log(pi1)
    offset <- pmax(a, b)
    denominator <- exp(a - offset) + exp(b - offset)
    list(r = exp(b - offset) / denominator, loglik = sum(offset + log(denominator)))
  }
  starts <- list()
  for (quantile in c(.05, .95)) for (initial_pi in c(.1, .4)) {
    mu <- unname(stats::quantile(x, quantile))
    tau2 <- 0
    pi1 <- initial_pi
    state <- posterior(mu, tau2, pi1)
    converged <- FALSE
    for (iteration in seq_len(opt$location_max_iter)) {
      r <- state$r
      pi1 <- min(1 - 1e-10, max(1e-10, mean(r)))
      mu <- sum(r * x / (variance + tau2)) / sum(r / (variance + tau2))
      objective <- function(v) sum(r * stats::dnorm(x, mu, sqrt(variance + v), log = TRUE))
      interior <- stats::optimize(objective, c(0, upper), maximum = TRUE, tol = 1e-7)$maximum
      candidates <- c(0, interior, upper)
      tau2 <- candidates[which.max(vapply(candidates, objective, numeric(1)))]
      next_state <- posterior(mu, tau2, pi1)
      increment <- next_state$loglik - state$loglik
      if (increment < -1e-7) stop("Location ECM likelihood decreased.")
      state <- next_state
      if (increment <= opt$location_tolerance) { converged <- TRUE; break }
    }
    starts[[length(starts) + 1L]] <- list(mean = mu * scale, variance = tau2 * scale^2,
      pi1 = pi1, loglik = state$loglik, converged = converged, iterations = iteration)
  }
  selected <- which.max(vapply(starts, function(x) x$loglik, numeric(1)))
  best <- starts[[selected]]
  best$enabled <- abs(best$mean) >= opt$minimum_location_se * scale
  best$selected_start <- selected
  best$start_loglik <- vapply(starts, function(x) x$loglik, numeric(1))
  best
}

adaptive_grid <- function(e, opt, locations) {
  r <- exp(seq(log(opt$r_min), log(opt$r_max), length.out = opt$J))
  make_axis <- function(se, location) {
    axis <- data.frame(mean = numeric(), variance = numeric(), j = integer(), type = character())
    if (opt$zero_anchor)
      axis <- rbind(axis, data.frame(mean = 0, variance = (r * se)^2, j = seq_along(r), type = "zero"))
    if (opt$learn_location && location$enabled) {
      rr <- if (opt$include_atom) c(0, r) else r
      jj <- if (opt$include_atom) c(0L, seq_along(r)) else seq_along(r)
      axis <- rbind(axis, data.frame(mean = location$mean, variance = (rr * se)^2, j = jj, type = "learned"))
    }
    if (!nrow(axis)) stop("No alternative basis remains; enable zero_anchor.")
    axis
  }
  a <- make_axis(median(sqrt(e$var_alpha)), locations$alpha)
  b <- make_axis(median(sqrt(e$var_beta)), locations$beta)
  component <- function(state, ma, mb, u, v, c, ja, jb, ta, tb, rho) {
    data.frame(state = state, mean_alpha = ma, mean_beta = mb,
      prior_var_alpha = u, prior_var_beta = v, prior_cov = c,
      scale_alpha = ja, scale_beta = jb, location_alpha = ta, location_beta = tb, rho = rho)
  }
  grid <- component("H00", 0, 0, 0, 0, 0, 0L, 0L, "null", "null", 0)
  grid <- rbind(grid, component("H10", a$mean, 0, a$variance, 0, 0, a$j, 0L, a$type, "null", 0),
                component("H01", 0, b$mean, 0, b$variance, 0, 0L, b$j, "null", b$type, 0))
  for (i in seq_len(nrow(a))) for (j in seq_len(nrow(b))) {
    if (opt$grid_mode == "same" && a$j[i] != b$j[j]) next
    rhos <- if (a$variance[i] == 0 || b$variance[j] == 0) 0 else opt$correlations
    for (rho in rhos)
      grid <- rbind(grid, component("H11", a$mean[i], b$mean[j], a$variance[i], b$variance[j],
                         rho * sqrt(a$variance[i] * b$variance[j]), a$j[i], b$j[j], a$type[i], b$type[j], rho))
  }
  grid$component_id <- paste0(grid$state, "_", seq_len(nrow(grid)))
  rownames(grid) <- NULL
  grid
}

adaptive_log_likelihood <- function(e, grid, covariance = ASH_BASE$ashmed_validate_estimates(e)) {
  logL <- matrix(0, nrow(e), nrow(grid))
  for (k in seq_len(nrow(grid))) {
    a <- e$var_alpha + grid$prior_var_alpha[k]
    b <- e$var_beta + grid$prior_var_beta[k]
    c <- covariance + grid$prior_cov[k]
    det <- a * b - c^2
    x <- e$alpha_hat - grid$mean_alpha[k]
    y <- e$beta_hat - grid$mean_beta[k]
    if (any(det <= 0 | !is.finite(det))) stop("Invalid adaptive marginal covariance.")
    logL[, k] <- -log(2 * pi) - .5 * (log(det) + (b * x^2 - 2 * c * x * y + a * y^2) / det)
  }
  stopifnot(all(is.finite(logL)))
  logL
}

adaptive_fit_weights <- function(logL, grid, opt) {
  if (opt$null_pseudocount == 0) return(ASH_BASE$ashmed_fit_weights(logL, grid, opt$base))
  offset <- apply(logL, 1L, max)
  F <- exp(logL - offset)
  p <- nrow(F)
  states <- c("H00", "H10", "H01", "H11")
  state_index <- match(grid$state, states)
  membership <- outer(state_index, seq_along(states), "==") * 1
  pseudo <- c(rep(opt$null_pseudocount / 3, 3), 0)
  total <- p + sum(pseudo)
  counts <- table(grid$state)
  w <- unname(opt$base$initial_pi[grid$state] / counts[grid$state])
  w <- as.numeric(w / sum(w))
  state_mass <- function(x) as.vector(crossprod(membership, x))
  objective <- function(x) {
    den <- as.vector(F %*% x)
    mass <- state_mass(x)
    if (any(den <= 0 | !is.finite(den)) || any(mass[1:3] <= 0)) return(-Inf)
    sum(log(den)) + sum(pseudo[1:3] * log(mass[1:3]))
  }
  gradient <- function(x) {
    mass <- state_mass(x)
    prior <- c(pseudo[1:3] / mass[1:3], 0)[state_index]
    (as.vector(crossprod(F, 1 / as.vector(F %*% x))) + prior) / total
  }
  em <- function(x) {
    ans <- x * gradient(x)
    ans / sum(ans)
  }
  value <- objective(w)
  trace <- data.frame(iteration = 0L, objective = value + sum(offset), dual_gap = NA_real_, relative_change = NA_real_)
  converged <- FALSE
  for (iteration in seq_len(opt$max_iter)) {
    w1 <- em(w)
    candidate <- w1
    if (opt$accelerate) {
      w2 <- em(w1)
      candidate <- w2
      r <- w1 - w
      v <- w2 - w1 - r
      if (sum(v^2) > 0 && sum(r^2) > 0) {
        step <- -max(1, min(1000, sqrt(sum(r^2) / sum(v^2))))
        base_value <- objective(w2)
        for (attempt in seq_len(15L)) {
          proposal <- pmax(0, w - 2 * step * r + step^2 * v)
          proposal <- proposal / sum(proposal)
          if (is.finite(objective(proposal))) {
            proposed <- em(proposal)
            if (objective(proposed) >= base_value) { candidate <- proposed; break }
          }
          step <- (step - 1) / 2
        }
      }
    }
    next_value <- objective(candidate)
    grad <- gradient(candidate)
    gap <- max(0, max(grad) - 1)
    if (gap > opt$kkt_tolerance &&
        ((next_value - value) / max(1, abs(value + sum(offset))) < opt$relative_tolerance || iteration %% 20L == 0L)) {
      best <- which.max(grad)
      vertex <- numeric(length(w)); vertex[best] <- 1
      line <- stats::optimize(function(t) -objective((1 - t) * candidate + t * vertex),
                              c(0, 1 - 1e-12), tol = 1e-10)
      if (-line$objective > next_value) {
        candidate <- (1 - line$minimum) * candidate + line$minimum * vertex
        next_value <- -line$objective
        gap <- max(0, max(gradient(candidate)) - 1)
      }
    }
    if (next_value < value - 1e-9 * max(1, abs(value))) stop("Adaptive penalized objective decreased.")
    relative <- abs(next_value - value) / max(1, abs(value + sum(offset)))
    w <- candidate / sum(candidate)
    value <- next_value
    trace[iteration + 1L, ] <- list(iteration, value + sum(offset), gap, relative)
    if (relative <= opt$relative_tolerance && gap <= opt$kkt_tolerance) { converged <- TRUE; break }
  }
  list(weights = w, loglik = sum(offset + log(as.vector(F %*% w))), objective = value + sum(offset),
       trace = trace, iterations = iteration, converged = converged, dual_gap = gap)
}

adaptive_moments <- function(e, grid, resp, covariance) {
  p <- nrow(e)
  ma <- mb <- sa <- sb <- ab <- aabb <- numeric(p)
  for (k in seq_len(nrow(grid))) {
    u <- grid$prior_var_alpha[k]; v <- grid$prior_var_beta[k]; c <- grid$prior_cov[k]
    a <- e$var_alpha + u; b <- e$var_beta + v; d <- covariance + c
    det <- a * b - d^2
    B11 <- (u * b - c * d) / det; B12 <- (-u * d + c * a) / det
    B21 <- (c * b - v * d) / det; B22 <- (-c * d + v * a) / det
    x <- e$alpha_hat - grid$mean_alpha[k]; y <- e$beta_hat - grid$mean_beta[k]
    am <- grid$mean_alpha[k] + B11 * x + B12 * y
    bm <- grid$mean_beta[k] + B21 * x + B22 * y
    A <- pmax(0, u - B11 * u - B12 * c); B <- pmax(0, v - B21 * c - B22 * v)
    C <- ((c - B11 * c - B12 * v) + (c - B21 * u - B22 * c)) / 2
    r <- resp[, k]
    ma <- ma + r * am; mb <- mb + r * bm
    sa <- sa + r * (A + am^2); sb <- sb + r * (B + bm^2)
    ab <- ab + r * (am * bm + C)
    aabb <- aabb + r * (am^2 * bm^2 + am^2 * B + bm^2 * A + 4 * am * bm * C + A * B + 2 * C^2)
  }
  data.frame(alpha_mean = ma, alpha_sd = sqrt(pmax(0, sa - ma^2)),
    beta_mean = mb, beta_sd = sqrt(pmax(0, sb - mb^2)), alpha_beta_cov = ab - ma * mb,
    mediation_mean = ab, mediation_sd = sqrt(pmax(0, aabb - ab^2)))
}

fit_adaptive_ashmed <- function(e, options = list()) {
  opt <- adaptive_options(options)
  covariance <- ASH_BASE$ashmed_validate_estimates(e)
  disabled <- list(mean = 0, enabled = FALSE, converged = TRUE, iterations = 0L)
  locations <- if (opt$learn_location) list(alpha = adaptive_location(e$alpha_hat, e$var_alpha, opt),
                                           beta = adaptive_location(e$beta_hat, e$var_beta, opt)) else
    list(alpha = disabled, beta = disabled)
  if (!opt$zero_anchor && !all(vapply(locations, function(x) x$enabled, logical(1)))) {
    # An unidentifiable near-zero location must not be labelled a nonnull atom.
    # Use the unchanged original shrinkage model when either learned mode
    # is insufficiently separated from zero. This rule uses estimates only.
    original <- ASH_BASE$fit_ashmed(e)
    grid <- original$components
    grid$mean_alpha <- grid$mean_beta <- 0
    return(list(posterior = original$posterior, lfdr = original$lfdr,
      moments = original$moments, grid = grid, pi = original$pi,
      locations = locations, fit = original$fit, options = opt, fallback = TRUE))
  }
  grid <- adaptive_grid(e, opt, locations)
  logL <- adaptive_log_likelihood(e, grid, covariance)
  fit <- adaptive_fit_weights(logL, grid, opt)
  resp <- ASH_BASE$ashmed_responsibilities(logL, fit$weights)
  states <- c("H00", "H10", "H01", "H11")
  posterior <- vapply(states, function(s) rowSums(resp[, grid$state == s, drop = FALSE]), numeric(nrow(e)))
  posterior <- matrix(posterior, nrow = nrow(e), dimnames = list(NULL, states))
  pi_state <- setNames(vapply(states, function(s) sum(fit$weights[grid$state == s]), numeric(1)), states)
  grid$weight <- fit$weights
  list(posterior = posterior, lfdr = pmin(1, pmax(0, rowSums(posterior[, 1:3, drop = FALSE]))),
       moments = adaptive_moments(e, grid, resp, covariance), grid = grid,
       pi = pi_state, locations = locations, fit = fit, options = opt, fallback = FALSE)
}

run_adaptive_ashmed <- function(input, q, options = list()) {
  fit <- fit_adaptive_ashmed(input$estimates, options)
  if (!fit$fit$converged) warning("Adaptive weight fit reached max_iter; inspect diagnostics.")
  if (!all(vapply(fit$locations, function(x) x$converged, logical(1))))
    warning("Adaptive marginal location fit reached max_iter; inspect diagnostics.")
  reject <- lfdr_step_up(fit$lfdr, q)
  ids <- if ("pathway_id" %in% names(input$estimates)) input$estimates$pathway_id else seq_along(reject)
  list(reject = reject, score = fit$lfdr, score_type = "local_FDR",
    posterior = data.frame(pathway_id = ids, fit$posterior, lfdr = fit$lfdr, reject = reject, fit$moments),
    components = fit$grid,
    diagnostics = list(engine = "adaptive_location_scale_ashmed", options = fit$options,
      locations = fit$locations, fallback = fit$fallback, pi = fit$pi, converged = fit$fit$converged,
      iterations = fit$fit$iterations, dual_gap = fit$fit$dual_gap,
      loglik = fit$fit$loglik, trace = fit$fit$trace,
      selected_mean_lfdr = if (any(reject)) mean(fit$lfdr[reject]) else NA_real_,
      likelihood = if (identical(input$meta$scenario, "binary")) "logistic_summary_normal_extension" else "linear_summary_normal"))
}

register_plugin <- function(registry) {
  # Cache fingerprint also records the private original-AS HMED dependency.
  register_method(registry, "ASHMED_Adaptive", run_adaptive_ashmed,
    version = paste0("1.0.0+base.", substr(file_hash(ASH_BASE_PATH), 1, 12)),
    description = "ASHMED extension: observed-data locations, independent scales, null-state Dirichlet penalty")
}
