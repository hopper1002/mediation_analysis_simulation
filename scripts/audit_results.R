# Read-only audit of the delivered demo; writes separate diagnostic artifacts.
# Run from project root: Rscript scripts/audit_results.R
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
audit_dir <- ensure_dir(file.path(PROJECT_ROOT, "docs", "audit"))
main_dir <- file.path(PROJECT_ROOT, "results", "run_45a379a50352")
source_dir <- file.path(PROJECT_ROOT, "results", "run_5580e7f5e7ad")
main_metrics <- utils::read.csv(file.path(main_dir, "replicate_metrics.csv"))
source_metrics <- utils::read.csv(file.path(source_dir, "replicate_metrics.csv"))
main_manifest <- read_manifest(file.path(PROJECT_ROOT, "data", "data_6fdb451c6d0d", "manifest.rds"))
source_manifest <- read_manifest(file.path(PROJECT_ROOT, "data", "data_9aaefd4678b4", "manifest.rds"))

saved_estimates <- function(data, row) {
  key <- object_hash(list(data_md5 = row$md5,
                          code = code_signature(PROJECT_ROOT, file.path(PROJECT_ROOT, "R", c("utils.R", "estimation.R"))),
                          R = R.version.string))
  path <- file.path(PROJECT_ROOT, "data", "estimates", paste0(key, ".rds"))
  if (file.exists(path)) readRDS(path)$estimates else estimate_coefficients(data)$estimates
}

oracle_lfdr <- function(data, e) {
  cfg <- data$dgp
  pi <- cfg$mixtures[[data$meta$mixture]]
  va <- if (data$meta$profile == "paper2026") cfg$alpha_noise_variance / data$meta$n else cfg$alpha_noise_variance^2
  vb <- if (data$meta$profile == "paper2026") cfg$beta_noise_variance / data$meta$n else cfg$beta_noise_variance^2
  a <- cfg$alpha_mean_scale * data$meta$tau
  b <- cfg$beta_mean_scale * data$meta$tau
  f0a <- stats::dnorm(e$alpha_hat, 0, sqrt(e$var_alpha))
  f1a <- stats::dnorm(e$alpha_hat, a, sqrt(e$var_alpha + va))
  f0b <- stats::dnorm(e$beta_hat, 0, sqrt(e$var_beta))
  f1b <- stats::dnorm(e$beta_hat, b, sqrt(e$var_beta + vb))
  # Uses population DGP parameters, solely to diagnose signal strength.
  # This Gaussian approximation with estimated SEs is not a calibrated oracle test.
  densities <- cbind(pi["H00"] * f0a * f0b, pi["H10"] * f1a * f0b,
                    pi["H01"] * f0a * f1b, pi["H11"] * f1a * f1b)
  rowSums(densities[, 1:3]) / rowSums(densities)
}

signal_rows <- list()
null_rows <- list()
selected_rows <- list()
datasets <- list()
for (which_run in c("main", "source")) {
  manifest <- if (which_run == "main") main_manifest else source_manifest
  metrics <- if (which_run == "main") main_metrics else source_metrics
  idx <- which(manifest$design$tau == 1.9 & manifest$design$n == 300)
  for (i in idx) {
    row <- manifest$design[i, ]
    data <- read_dataset(manifest, row$case_id)
    e <- saved_estimates(data, row)
    truth <- data$truth
    alpha_active <- truth$state %in% c("H10", "H11")
    h11 <- truth$is_nonnull
    lf <- oracle_lfdr(data, e)
    signal_rows[[length(signal_rows) + 1L]] <- data.frame(
      run = which_run, row[, c("case_id", "scenario", "mixture", "n", "tau", "replicate")],
      X_ones = sum(data$X), alpha_mean = data$dgp$alpha_mean_scale * row$tau,
      alpha_prior_sd = if (row$profile == "paper2026") sqrt(data$dgp$alpha_noise_variance / row$n) else data$dgp$alpha_noise_variance,
      median_alpha_SE = median(e$se_alpha),
      mean_abs_true_alpha_over_SE = mean(abs(truth$alpha[alpha_active]) / e$se_alpha[alpha_active]),
      alpha_p05_among_H11 = mean(e$p_alpha[h11] < 0.05),
      beta_p05_among_H11 = mean(e$p_beta[h11] < 0.05),
      both_p05_among_H11 = mean(e$p_alpha[h11] < 0.05 & e$p_beta[h11] < 0.05),
      approximate_oracle_min_lfdr = min(lf),
      approximate_oracle_discoveries = sum(lfdr_step_up(lf, 0.05)))
    if (which_run == "main" && row$scenario == "binary" && row$mixture == "sparse") {
      hdmt_row <- metrics[metrics$case_id == row$case_id & metrics$method == "HDMT", ]
      outcome <- readRDS(file.path(PROJECT_ROOT, "results", "method_cache", "HDMT", paste0(hdmt_row$cache_key, ".rds")))$result
      fitted_null <- outcome$diagnostics$null
      true_pi <- table(factor(truth$state, levels = c("H00", "H10", "H01", "H11"))) / nrow(truth)
      true_null <- list(alpha00 = unname(true_pi["H00"]), alpha01 = unname(true_pi["H01"]),
                        alpha10 = unname(true_pi["H10"]), alpha1 = unname(true_pi["H00"] + true_pi["H01"]),
                        alpha2 = unname(true_pi["H00"] + true_pi["H10"]))
      pv <- as.matrix(e[, c("p_alpha", "p_beta")])
      fdr_oracle <- do.call(HDMT::fdr_est, c(true_null, list(input_pvalues = pv, exact = 0)))
      p_max <- pmax(e$p_alpha, e$p_beta)
      eligible <- fdr_oracle <= 0.05
      rej_oracle <- if (any(eligible)) p_max <= max(p_max[eligible]) else rep(FALSE, nrow(e))
      null_rows[[length(null_rows) + 1L]] <- data.frame(
        case_id = row$case_id, replicate = row$replicate,
        true_H00 = true_pi["H00"], true_H10 = true_pi["H10"], true_H01 = true_pi["H01"], true_H11 = true_pi["H11"],
        estimated_H00 = fitted_null$alpha00, estimated_H10 = fitted_null$alpha10,
        estimated_H01 = fitted_null$alpha01,
        estimated_H11 = 1 - fitted_null$alpha00 - fitted_null$alpha10 - fitted_null$alpha01,
        discoveries = sum(outcome$reject), true_null_proportions_diagnostic_R = sum(rej_oracle))
      if (any(outcome$reject)) {
        sel <- which(outcome$reject)
        selected_rows[[length(selected_rows) + 1L]] <- data.frame(
          case_id = row$case_id, replicate = row$replicate,
          truth[sel, c("pathway_id", "state", "alpha", "beta")],
          e[sel, c("alpha_hat", "beta_hat", "se_alpha", "se_beta", "p_alpha", "p_beta")],
          estimated_FDR = outcome$score[sel], true_null_proportions_estimated_FDR = fdr_oracle[sel])
      }
    }
    if (row$replicate == 1 && (row$scenario == "linear" || (which_run == "main" && row$scenario == "binary" && row$mixture == "sparse"))) {
      datasets[[paste(which_run, row$case_id)]] <- list(data = data, e = e)
    }
  }
}
signals <- do.call(rbind, signal_rows)
nulls <- do.call(rbind, null_rows)
selected <- do.call(rbind, selected_rows)
write_csv(signals, file.path(audit_dir, "signal_strength.csv"))
write_csv(nulls, file.path(audit_dir, "null_proportions.csv"))
write_csv(selected, file.path(audit_dir, "false_positive_paths.csv"))
problem <- main_metrics[main_metrics$scenario == "binary" & main_metrics$mixture == "sparse" & main_metrics$n == 300 & main_metrics$tau == 1.9,
                        c("case_id", "replicate", "method", "discoveries", "false_discoveries", "true_discoveries", "FDP", "power")]
write_csv(problem, file.path(audit_dir, "problem_condition.csv"))
cat("Signal audit by run/scenario/mixture:\n")
print(stats::aggregate(cbind(median_alpha_SE, mean_abs_true_alpha_over_SE, alpha_p05_among_H11, beta_p05_among_H11, approximate_oracle_discoveries) ~ run + scenario + mixture, signals, mean), row.names = FALSE)
cat("Null proportion diagnosis:\n")
print(nulls, row.names = FALSE)
cat("False positive pathways:\n")
print(selected, row.names = FALSE)

# Isolate optimization bounds without changing initialization or convergence tolerance.
em_rows <- list()
em_fit <- getFromNamespace("EM_fun", "MLFDR")
for (name in names(datasets)) {
  ds <- datasets[[name]]
  e <- ds$e
  d <- ds$data
  for (bounds in c("package_default", "wide_bounds")) {
    cat("EM audit:", name, bounds, "\n")
    default_kappa <- c(min(0.1, min(e$var_alpha)), max(10, max(e$var_alpha)))
    default_psi <- c(min(0.1, min(e$var_beta)), max(10, max(e$var_beta)))
    fit <- em_fit(cbind(e$alpha_hat, e$beta_hat), var_alpha = e$var_alpha, var_beta = e$var_beta,
                  kappa.init = 1, psi.init = 1,
                  kappa_int = if (bounds == "wide_bounds") c(1e-8, 100) else NULL,
                  psi_int = if (bounds == "wide_bounds") c(1e-8, 100) else NULL,
                  epsilon = 0.01, verbose = FALSE)
    densities <- vapply(seq_len(4), function(j) fit$lambda[j] * stats::dnorm(e$alpha_hat, fit$mu[[j]][1], sqrt(fit$sigma[j, , 1, 1])) *
                          stats::dnorm(e$beta_hat, fit$mu[[j]][2], sqrt(fit$sigma[j, , 2, 2])), numeric(nrow(e)))
    lfdr <- rowSums(densities[, 1:3]) / rowSums(densities)
    cm <- compute_metrics(lfdr_step_up(lfdr, 0.05), d$truth)
    em_rows[[length(em_rows) + 1L]] <- data.frame(
      run = strsplit(name, " ", fixed = TRUE)[[1]][1], case_id = d$meta$case_id, bounds = bounds,
      true_alpha_variance = if (d$meta$profile == "paper2026") 1 / d$meta$n else 1,
      true_beta_variance = if (d$meta$profile == "paper2026") 4 / d$meta$n else 16,
      kappa_lower = if (bounds == "wide_bounds") 1e-8 else default_kappa[1],
      psi_upper = if (bounds == "wide_bounds") 100 else default_psi[2],
      fitted_kappa = fit$sigma[2, 1, 1, 1] - e$var_alpha[1],
      fitted_psi = fit$sigma[3, 1, 2, 2] - e$var_beta[1],
      fitted_mu = fit$mu[[4]][1], fitted_theta = fit$mu[[4]][2], fitted_pi11 = fit$lambda[4],
      loglik = fit$loglik, as.data.frame(cm))
    write_csv(do.call(rbind, em_rows), file.path(audit_dir, "em_bounds_diagnostic.csv"))
  }
}
cat("EM bounds diagnostics (single datasets, not a new comparison experiment):\n")
print(do.call(rbind, em_rows), row.names = FALSE)
cat("Artifacts:", audit_dir, "\n")
