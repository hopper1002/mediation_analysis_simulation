# Validate all four methods against all saved raw datasets and caches.
# No generation and no modification of experiment outputs.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source(file.path(PROJECT_ROOT, "R/bootstrap.R"), encoding = "UTF-8")
source(file.path(PROJECT_ROOT, "R/ashmed_comparison.R"), encoding = "UTF-8")
pointer_path <- file.path(PROJECT_ROOT, "docs/current_ashmed.rds")
if (!file.exists(pointer_path)) stop("Run notebook 03 or scripts/run_ashmed.R first.")
pointer <- readRDS(pointer_path)
inputs <- read_ashmed_inputs(PROJECT_ROOT)
eq <- function(x, y, tolerance = 1e-11) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
}
compare_tables <- function(x, y) {
  stopifnot(identical(names(x), names(y)), nrow(x) == nrow(y))
  for (name in names(x)) eq(x[[name]], y[[name]])
}
read_csv_like <- function(path, reference) {
  classes <- vapply(reference, function(x) if (is.character(x)) "character" else
    if (is.logical(x)) "logical" else "numeric", character(1))
  utils::read.csv(path, stringsAsFactors = FALSE, colClasses = classes)
}
report <- character()
total_cases <- total_calls <- 0L
for (label in names(pointer)) {
  dir <- pointer[[label]]$run_dir
  config <- readRDS(file.path(dir, "config.rds"))
  manifest <- read_manifest(pointer[[label]]$manifest_path)
  saved <- readRDS(file.path(dir, "experiment.rds"))
  experiment <- list(metrics = saved$metrics, summary = saved$summary,
                     config = config, manifest = manifest, output_dir = dir)
  stopifnot(identical(manifest$dgp, dgp_config(config)),
            identical(config$methods, c("HDMT", "MDACT", "MLFDR", "ASHMED")),
            nrow(manifest$design) == 720L, nrow(experiment$metrics) == 2880L,
            all(experiment$metrics$status == "ok"), config$q == .05,
            all(config$n == 300L), config$m == 1000L, config$repetitions == 40L)
  csv_metrics <- read_csv_like(file.path(dir, "replicate_metrics.csv"), experiment$metrics)
  compare_tables(experiment$metrics, csv_metrics)
  summary <- summarise_metrics(experiment$metrics)
  paired <- summarise_paired_metrics(experiment$metrics, left = "ASHMED", right = c("HDMT", "MDACT", "MLFDR"))
  compare_tables(summary, read_csv_like(file.path(dir, "summary.csv"), summary))
  compare_tables(paired, read_csv_like(file.path(dir, "paired_differences.csv"), paired))
  audit <- audit_ashmed_baseline(experiment, inputs[[label]])
  stopifnot(all(audit$unchanged))
  recorded_data_audit <- utils::read.csv(file.path(dir, "data_audit.csv"))
  stopifnot(recorded_data_audit$checksum_matches == 720L,
            recorded_data_audit$mtime_unchanged, recorded_data_audit$manifest_unchanged)
  fits <- collect_ashmed_fit_diagnostics(experiment, PROJECT_ROOT)
  compare_tables(fits, read_csv_like(file.path(dir, "ashmed_fit_diagnostics.csv"), fits))
  stopifnot(nrow(fits) == 720L, all(fits$converged), all(fits$n_fit == 1000L),
            all(fits$n_components == 61L), all(fits$dual_gap <= 1e-4 + 1e-12))
  for (i in seq_len(nrow(manifest$design))) {
    case_id <- manifest$design$case_id[i]
    data <- read_dataset(manifest, case_id) # verifies the original raw MD5
    rows <- experiment$metrics[experiment$metrics$case_id == case_id, ]
    results <- get_case_results(experiment, case_id, PROJECT_ROOT)
    for (method in config$methods) {
      result <- validate_method_result(results[[method]], ncol(data$M))
      actual <- compute_metrics(result$reject, data$truth)
      row <- rows[rows$method == method, ]
      for (name in names(actual)) eq(actual[[name]], row[[name]])
    }
    result <- results$ASHMED
    posterior <- result$posterior
    stopifnot(identical(posterior$pathway_id, colnames(data$M)),
              identical(result$reject, posterior$reject),
              all(is.finite(as.matrix(posterior[, setdiff(names(posterior), "pathway_id")]))))
    probabilities <- as.matrix(posterior[, c("H00", "H10", "H01", "H11")])
    stopifnot(all(probabilities >= -1e-12 & probabilities <= 1 + 1e-12),
              max(abs(rowSums(probabilities) - 1)) < 1e-12,
              all(posterior$alpha_sd >= 0 & posterior$beta_sd >= 0 & posterior$mediation_sd >= 0))
    eq(posterior$lfdr, result$score)
    eq(posterior$lfdr, rowSums(probabilities[, 1:3]))
    eq(posterior$lfdr, 1 - probabilities[, 4])
    stopifnot(identical(lfdr_step_up(posterior$lfdr, config$q), result$reject),
              !any(result$reject) || mean(posterior$lfdr[result$reject]) <= config$q + 1e-12,
              all(diff(result$diagnostics$trace$loglik) >= -1e-8),
              result$diagnostics$converged)
    components <- result$components
    eq(sum(components$weight), 1)
    for (state in c("H00", "H10", "H01", "H11")) {
      subset <- components$state == state
      eq(sum(components$weight[subset]), result$diagnostics$pi[state])
      eq(sum(components$within_state_weight[subset]), 1)
    }
    stopifnot(all(components$weight >= 0),
              all(components$prior_var_beta[components$state == "H10"] == 0),
              all(components$prior_var_alpha[components$state == "H01"] == 0))
    expected <- if (manifest$design$scenario[i] == "binary") "logistic_summary_normal_extension" else "linear_summary_normal"
    stopifnot(result$diagnostics$likelihood == expected)
    if (i %% 120L == 0L) message(label, ": validated ", i, "/720 saved datasets")
  }
  total_cases <- total_cases + nrow(manifest$design)
  total_calls <- total_calls + nrow(experiment$metrics)
  report <- c(report, sprintf("%s: 720 saved datasets, 2880 successful calls; ASHMED 720/720 converged; max gap %.6g; old-method differences <= %.3g.",
                              label, max(fits$dual_gap), max(audit$max_absolute_difference)))
}
report <- c("ASHMED delivery validation passed.", report,
            sprintf("Total: %d saved datasets; %d method results checked against truth and rejection caches.", total_cases, total_calls),
            "Original data checksums match; baseline R/V/S/A/FDP/power preserved; summary/paired tables recomputed and matched.",
            "No experiment dataset generation was invoked.")
writeLines(report, file.path(PROJECT_ROOT, "docs/ASHMED_VALIDATION.txt"), useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n")
