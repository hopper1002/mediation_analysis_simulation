run_hdmt <- function(input, q, options = list()) {
  e <- input$estimates
  p <- as.matrix(e[, c("p_alpha", "p_beta")])
  null <- HDMT::null_estimation(p)
  fdr <- HDMT::fdr_est(null$alpha00, null$alpha01, null$alpha10,
                       null$alpha1, null$alpha2, p, exact = if (is.null(options$exact)) 0L else options$exact)
  if (length(fdr) != nrow(e) || any(!is.finite(fdr))) stop("HDMT returned invalid estimated FDR.")
  p_max <- pmax(e$p_alpha, e$p_beta)
  eligible <- which(fdr <= q)
  threshold <- if (length(eligible)) max(p_max[eligible]) else -Inf
  list(reject = p_max <= threshold, score = fdr, score_type = "estimated_FDR",
       diagnostics = list(null = null, threshold = threshold))
}

run_mdact <- function(input, q, options = list()) {
  e <- input$estimates
  null <- HDMT::null_estimation(as.matrix(e[, c("p_alpha", "p_beta")]))
  idx <- mdact_vendor$balancing_DACT_control_DR_adjust(
    e$p_alpha, e$p_beta, null, significance_upper = q, control.method = "FDR")
  if (anyNA(idx) || any(idx != as.integer(idx)) || any(idx < 1L | idx > nrow(e)))
    stop("MDACT returned invalid rejection indices.")
  reject <- rep(FALSE, nrow(e))
  reject[idx] <- TRUE
  list(reject = reject, diagnostics = list(null = null,
       implementation = "vendor/mdact_core.R; original threshold search, equivalent analytic F00"))
}

bound_lfdr <- function(lfdr, tolerance = 1e-12) {
  stopifnot(is.numeric(lfdr), length(lfdr) > 0L, all(is.finite(lfdr)))
  if (any(lfdr < -tolerance | lfdr > 1 + tolerance))
    stop("MLFDR returned probabilities outside [0,1] beyond numeric tolerance.")
  pmin(1, pmax(0, lfdr))
}

# Replaces the package's k=1 / while(k<m) edge cases (zero or all rejections).
# Cutoffs are considered at tied-score group ends, so ties cannot inflate the selected average.
lfdr_step_up <- function(lfdr, q) {
  lfdr <- bound_lfdr(lfdr)
  stopifnot(q > 0, q < 1)
  ord <- order(lfdr)
  values <- lfdr[ord]
  group_ends <- c(which(diff(values) != 0), length(values))
  averages <- cumsum(values) / seq_along(values)
  ok <- group_ends[averages[group_ends] <= q]
  reject <- rep(FALSE, length(values))
  if (length(ok)) reject[ord[seq_len(max(ok))]] <- TRUE
  reject
}

run_mlfdr <- function(input, q, options = list()) {
  e <- input$estimates
  engine <- if (is.null(options$engine)) "package" else options$engine
  if (engine == "paper_em") {
    options$engine <- NULL
    fit <- fit_mlfdr_paper(e$alpha_hat, e$beta_hat, e$var_alpha, e$var_beta, input$meta$n, options)
    lfdr <- bound_lfdr(fit$lfdr)
    reject <- lfdr_step_up(lfdr, q)
    return(list(reject = reject, score = lfdr, score_type = "local_FDR",
                diagnostics = list(engine = "paper_em", pi = setNames(fit$pi, c("H00", "H10", "H01", "H11")),
                                   mu = fit$mu, theta = fit$theta, kappa = fit$kappa, psi = fit$psi,
                                   loglik = fit$loglik, iterations = fit$iterations, converged = fit$converged,
                                   start_loglik = fit$start_loglik, start_converged = fit$start_converged,
                                   selected_start = fit$selected_start, scale = fit$scale,
                                   variance_upper = fit$upper,
                                   roundoff_clipped = sum(lfdr != fit$lfdr),
                                   selected_mean_lfdr = if (any(reject)) mean(lfdr[reject]) else NA_real_)))
  }
  if (engine != "package") stop("Unknown MLFDR engine: ", engine)
  options$engine <- NULL
  defaults <- list(eps = 0.01, twostep = FALSE, verbose = FALSE)
  opts <- utils::modifyList(defaults, options)
  forbidden <- intersect(names(opts), c("alpha", "beta", "var_alpha", "var_beta"))
  if (length(forbidden)) stop("MLFDR options cannot override input coefficients.")
  fit <- do.call(MLFDR::localFDR, c(list(alpha = e$alpha_hat, beta = e$beta_hat,
                                       var_alpha = e$var_alpha, var_beta = e$var_beta), opts))
  raw_lfdr <- fit$lfdr
  lfdr <- bound_lfdr(raw_lfdr)
  reject <- lfdr_step_up(lfdr, q)
  list(reject = reject, score = lfdr, score_type = "local_FDR",
       diagnostics = list(pi = fit$pi, roundoff_clipped = sum(lfdr != raw_lfdr),
                          selected_mean_lfdr = if (any(reject)) mean(lfdr[reject]) else NA_real_))
}
