# Check one existing replicate in every condition of both saved datasets.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source(file.path(PROJECT_ROOT, "R/bootstrap.R"), encoding = "UTF-8")
registry <- load_method_registry(PROJECT_ROOT)
sources <- list(paper_adjusted = c("data_f87521c2b4eb", "run_2836e6c08e18"),
                teaching_fixed = c("data_c054ed713a76", "run_af52854e0ae0"))
output <- list()
for (label in names(sources)) {
  spec <- sources[[label]]
  manifest <- read_manifest(file.path(PROJECT_ROOT, "data", spec[1], "manifest.rds"))
  config <- readRDS(file.path(PROJECT_ROOT, "results", spec[2], "config.rds"))
  config$methods <- c("HDMT", "MDACT", "MLFDR", "ASHMED")
  config$method_options$ASHMED <- list()
  idx <- which(manifest$design$replicate == 1L)
  rows <- lapply(idx, function(i) analyse_case(i, manifest, config, PROJECT_ROOT, registry))
  metrics <- do.call(rbind, rows)
  metrics$dataset <- label
  summary <- summarise_metrics(metrics)
  summary$dataset <- label
  fits <- lapply(seq_len(nrow(metrics)), function(i) {
    row <- metrics[i, ]
    if (row$method != "ASHMED" || row$status != "ok") return(NULL)
    result <- readRDS(file.path(PROJECT_ROOT, "results/method_cache/ASHMED", paste0(row$cache_key, ".rds")))$result
    data.frame(case_id = row$case_id, dataset = label, converged = result$diagnostics$converged,
               iterations = result$diagnostics$iterations, dual_gap = result$diagnostics$dual_gap)
  })
  write_csv(metrics, file.path(PROJECT_ROOT, "tmp", paste0("ashmed_pilot_", label, ".csv")))
  print(summary[, c("dataset", "scenario", "mixture", "tau", "method", "FDR_mean", "power_mean")], row.names = FALSE)
  print(do.call(rbind, fits), row.names = FALSE)
  stopifnot(all(metrics$status == "ok"), all(do.call(rbind, fits)$converged),
            all(metrics$cached[metrics$method != "ASHMED"]))
}
