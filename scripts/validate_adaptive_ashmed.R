# Independently recompute final metrics from saved truth and all five decision caches.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source("R/bootstrap.R", encoding = "UTF-8")
source("R/ashmed_optimization.R", encoding = "UTF-8")
campaign <- optimization_directory(PROJECT_ROOT)
lock <- readRDS(file.path(campaign, "locked_method.rds"))
stopifnot(identical(method_fingerprint(load_method_registry(PROJECT_ROOT)$ASHMED_Adaptive, PROJECT_ROOT),
                    lock$selection$fingerprint), all(lock$confirmation$passes),
          !length(intersect(lock$selection$development_replicates, lock$selection$validation_replicates)),
          !length(intersect(lock$selection$confirmation_replicates, lock$selection$validation_replicates)))
protected <- audit_protected_adaptive_files(PROJECT_ROOT)
inputs <- optimization_inputs(PROJECT_ROOT)
pointer <- readRDS(file.path(PROJECT_ROOT, "docs/current_ashmed_adaptive.rds"))
eq <- function(x, y, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}
compare_csv <- function(expected, path) {
  classes <- vapply(expected, function(x) if (is.character(x)) "character" else
    if (is.logical(x)) "logical" else "numeric", character(1))
  actual <- read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
  stopifnot(identical(names(expected), names(actual)), nrow(expected) == nrow(actual))
  for (name in names(expected)) eq(expected[[name]], actual[[name]])
}
report <- "ASHMED_Adaptive final evaluation validation passed."
total_cases <- total_calls <- 0L
for (label in names(pointer)) {
  dir <- pointer[[label]]$directory
  x <- readRDS(file.path(dir, "comparison.rds"))
  e <- x$experiment
  manifest <- e$manifest
  metrics <- x$combined
  stopifnot(nrow(manifest$design) == 504L, nrow(metrics) == 2520L,
    identical(sort(unique(manifest$design$replicate)), 13:40),
    identical(manifest$dgp, inputs[[label]]$manifest$dgp),
    identical(manifest$directory, inputs[[label]]$manifest$directory),
    identical(e$config$methods, "ASHMED_Adaptive"), e$config$q == .05,
    all(metrics$status == "ok"), all(is.finite(metrics$FDP)), all(is.finite(metrics$power)),
    all(metrics$FDP >= 0 & metrics$FDP <= 1), all(metrics$power >= 0 & metrics$power <= 1))
  canonical <- environment(load_method_registry(PROJECT_ROOT)$ASHMED_Adaptive$run)$adaptive_options(lock$selection$options)
  canonical$base <- NULL
  stopifnot(identical(e$config$method_options$ASHMED_Adaptive, canonical))
  compare_csv(metrics, file.path(dir, "replicate_metrics.csv"))
  s <- summarise_metrics(metrics); s$dataset <- label; s$stage <- "validation"; s$arm <- lock$selection$arm
  compare_csv(s, file.path(dir, "summary.csv"))
  compare_csv(adaptive_macro_summary(metrics), file.path(dir, "macro_summary.csv"))
  paired <- summarise_paired_metrics(metrics, "ASHMED_Adaptive", c("HDMT", "MDACT", "MLFDR", "ASHMED"))
  compare_csv(paired, file.path(dir, "paired_differences.csv"))
  fits <- collect_adaptive_diagnostics(e, PROJECT_ROOT)
  compare_csv(fits, file.path(dir, "adaptive_fit_diagnostics.csv"))
  stopifnot(all(fits$converged), all(fits$location_converged), all(fits$dual_gap <= 1e-4 + 1e-12))
  old <- inputs[[label]]$reference
  old <- old[old$case_id %in% manifest$design$case_id, ]
  new <- metrics[metrics$method != "ASHMED_Adaptive", ]
  key <- function(d) paste(d$case_id, d$method, sep = "::")
  stopifnot(setequal(key(old), key(new)), !anyDuplicated(key(new)), all(new$cached))
  new <- new[match(key(old), key(new)), ]
  for (col in c("discoveries", "false_discoveries", "true_discoveries", "alternatives", "FDP", "power", "cache_key"))
    eq(new[[col]], old[[col]])
  # Also verify all 720 raw files, including development repeats, against original manifest MD5.
  original <- inputs[[label]]$manifest
  stopifnot(all(unname(tools::md5sum(file.path(original$directory, original$design$file))) == original$design$md5))
  for (i in seq_len(nrow(manifest$design))) {
    id <- manifest$design$case_id[i]
    data <- read_dataset(manifest, id)
    rows <- metrics[metrics$case_id == id, ]
    stopifnot(nrow(rows) == 5L)
    for (j in seq_len(nrow(rows))) {
      row <- rows[j, ]
      outcome <- readRDS(file.path(PROJECT_ROOT, "results/method_cache", row$method, paste0(row$cache_key, ".rds")))
      result <- validate_method_result(outcome$result, ncol(data$M))
      actual <- compute_metrics(result$reject, data$truth)
      for (name in names(actual)) eq(actual[[name]], row[[name]])
      if (row$method != "ASHMED_Adaptive") next
      p <- result$posterior; d <- result$diagnostics; grid <- result$components
      probabilities <- as.matrix(p[, c("H00", "H10", "H01", "H11")])
      stopifnot(identical(p$pathway_id, colnames(data$M)), identical(p$reject, result$reject),
        all(is.finite(as.matrix(p[, setdiff(names(p), "pathway_id")]))),
        all(probabilities >= -1e-12 & probabilities <= 1 + 1e-12),
        max(abs(rowSums(probabilities) - 1)) < 1e-12,
        all(p$alpha_sd >= 0 & p$beta_sd >= 0 & p$mediation_sd >= 0),
        identical(lfdr_step_up(p$lfdr, .05), result$reject),
        !any(result$reject) || mean(p$lfdr[result$reject]) <= .05 + 1e-12,
        d$converged, all(vapply(d$locations, function(z) z$converged, logical(1))))
      eq(p$lfdr, result$score); eq(p$lfdr, rowSums(probabilities[, 1:3])); eq(p$lfdr, 1 - probabilities[, 4])
      trace <- d$trace
      target <- if ("objective" %in% names(trace)) trace$objective else trace$loglik
      stopifnot(all(diff(target) >= -1e-7), all(grid$weight >= 0))
      eq(sum(grid$weight), 1)
      for (state in names(d$pi)) eq(sum(grid$weight[grid$state == state]), d$pi[state])
      stopifnot(all(grid$prior_var_beta[grid$state == "H10"] == 0),
                all(grid$prior_var_alpha[grid$state == "H01"] == 0),
                nrow(grid) == if (d$fallback) 61L else 64L)
      expected <- if (manifest$design$scenario[i] == "binary") "logistic_summary_normal_extension" else "linear_summary_normal"
      stopifnot(d$likelihood == expected)
    }
    if (i %% 100L == 0L) message(label, ": ", i, "/504 data validated")
  }
  total_cases <- total_cases + 504L; total_calls <- total_calls + 2520L
  report <- c(report, sprintf("%s: 504 held-out cases / 2520 five-method results; all adaptive fits converged; max gap %.6g; fallbacks %d; raw MD5 720/720.",
    label, max(fits$dual_gap), sum(fits$fallback)))
}
report <- c(report, sprintf("Total: %d saved cases and %d method results recomputed from truth and decision caches.", total_cases, total_calls),
  "All summaries, MCSE/interval and paired tables recomputed and matched CSV.",
  "Frozen method fingerprint and options match development selection; repeats 13:40 disjoint from development/confirmation.",
  sprintf("Protected original files: %d/%d unchanged.", sum(protected$unchanged), nrow(protected)),
  "No original simulation data regeneration; baseline decisions and metrics exactly preserved.")
writeLines(report, file.path(PROJECT_ROOT, "docs/ASHMED_ADAPTIVE_VALIDATION.txt"), useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n")
