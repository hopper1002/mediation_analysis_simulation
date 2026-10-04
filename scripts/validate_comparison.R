PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
pointer <- readRDS(file.path(PROJECT_ROOT, "docs", "current_comparison.rds"))
manifest <- read_manifest(pointer$manifest_path)
run_dir <- pointer$output_dir
config <- pointer$config
rebuilt <- build_design(config)
stopifnot(identical(rebuilt$case_id, manifest$design$case_id),
          identical(rebuilt$seed, manifest$design$seed))
generator <- code_signature(PROJECT_ROOT, file.path(PROJECT_ROOT, "R", c("utils.R", "config.R", "simulation.R")))
expected_data_id <- paste0("data_", substr(object_hash(list(dgp=dgp_config(config), generator=generator)), 1, 12))
stopifnot(identical(expected_data_id, manifest$data_id))
record <- readRDS(file.path(run_dir, "experiment.rds"))
metrics <- record$metrics
stopifnot(nrow(metrics) == nrow(manifest$design) * length(config$methods),
          all(metrics$status == "ok"), all(metrics$q == config$q),
          all(metrics$FDP >= 0 & metrics$FDP <= 1), all(metrics$power >= 0 & metrics$power <= 1),
          !anyDuplicated(paste(metrics$case_id, metrics$method)),
          identical(manifest$dgp, dgp_config(config)))
for (i in seq_len(nrow(manifest$design))) {
  id <- manifest$design$case_id[i]
  data <- read_dataset(manifest, id)
  rows <- metrics[metrics$case_id == id, ]
  for (j in seq_len(nrow(rows))) {
    row <- rows[j, ]
    path <- file.path(PROJECT_ROOT, "results", "method_cache", row$method, paste0(row$cache_key, ".rds"))
    outcome <- readRDS(path)
    actual <- compute_metrics(outcome$result$reject, data$truth)
    for (name in names(actual)) stopifnot(isTRUE(all.equal(as.numeric(row[[name]]), as.numeric(actual[[name]]), tolerance=1e-12)))
  }
}
stopifnot(isTRUE(all.equal(summarise_metrics(metrics), record$summary, tolerance=1e-12)))
paired <- utils::read.csv(file.path(run_dir, "paired_differences.csv"))
stopifnot(isTRUE(all.equal(summarise_paired_metrics(metrics), paired, tolerance=1e-12, check.attributes=FALSE)))
fits <- utils::read.csv(file.path(run_dir, "fit_diagnostics.csv"))
stopifnot(nrow(fits) == nrow(manifest$design), all(fits$converged))
provenance <- readRDS(file.path(run_dir, "provenance.rds"))
stopifnot(identical(provenance$code, code_signature(PROJECT_ROOT)))
for (scenario in config$scenarios) stopifnot(file.exists(file.path(run_dir, "figures", paste0(scenario, "_comparison.png"))))
# This validator runs through Rscript after IRkernel generated the notebook:
# identical case seeds/data id and a cache-only call establish entry-point parity.
probe <- analyse_case(1L, manifest, config, PROJECT_ROOT, load_method_registry(PROJECT_ROOT))
stopifnot(all(probe$cached), all(probe$status == "ok"))
writeLines(c(paste0("run_dir=", run_dir), paste0("manifest_path=", pointer$manifest_path)),
            file.path(PROJECT_ROOT, "docs", "current_comparison.txt"))
cat(sprintf("PASS: %d saved datasets and %d method decisions/metrics; summaries, paired differences, all EM convergence, source provenance, figures, Rscript/IRkernel seed and cache parity.\n",
            nrow(manifest$design), nrow(metrics)))
cat("Run:", run_dir, "\n")
