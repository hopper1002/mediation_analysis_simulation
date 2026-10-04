# Small end-to-end integration check, independent from the delivered demo.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash = "/")
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
# Keep integration-test datasets and deliberately failed runs out of production results.
fixture_root <- ensure_dir(file.path(PROJECT_ROOT, "tmp", "pipeline_fixture"))
for (folder in c("R", "vendor")) {
  from <- list.files(file.path(PROJECT_ROOT, folder), recursive = TRUE, full.names = TRUE)
  for (path in from) {
    target <- file.path(fixture_root, substring(path, nchar(PROJECT_ROOT) + 2L))
    ensure_dir(dirname(target))
    if (!file.copy(path, target, overwrite = TRUE)) stop("Cannot prepare test fixture: ", path)
  }
}
PROJECT_ROOT <- fixture_root
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
config <- default_config()
config$m <- 200L
config$n <- 80L
config$tau <- 1.9
config$scenarios <- "linear"
config$mixtures <- config$mixtures["dense"]
config$repetitions <- 2L
config$methods <- c("HDMT", "MDACT", "MLFDR")
config$export_first_long_csv <- FALSE
config$workers <- 1L
manifest <- generate_data_files(config, PROJECT_ROOT, verbose = FALSE)
manifest <- read_manifest(file.path(manifest$directory, "manifest.rds"))
one <- read_dataset(manifest, manifest$design$case_id[1])
stopifnot(identical(one, readRDS(file.path(manifest$directory, manifest$design$file[1]))))
serial <- run_experiment(manifest, config, PROJECT_ROOT)
stopifnot(nrow(serial$metrics) == 6L, all(serial$metrics$status == "ok"))
config$workers <- 2L
parallel <- run_experiment(manifest, config, PROJECT_ROOT)
columns <- c("case_id", "method", "FDP", "power", "discoveries")
stopifnot(identical(serial$metrics[, columns], parallel$metrics[, columns]), all(parallel$metrics$cached))
cat("PASS: file round-trip, serial/PSOCK equivalence, successful-result cache reuse.\n")

config$methods <- c(config$methods, "MaxP_BH")
extended <- run_experiment(manifest, config, PROJECT_ROOT)
stopifnot(nrow(extended$metrics) == 8L, all(extended$metrics$status == "ok"),
          all(extended$metrics$cached[extended$metrics$method != "MaxP_BH"]))
cat("PASS: adding a plugin reuses existing data and unchanged method results.\n")

config$q <- 0.06
changed_q <- run_experiment(manifest, config, PROJECT_ROOT)
stopifnot(all(changed_q$metrics$status == "ok"),
          !any(changed_q$metrics$cache_key %in% extended$metrics$cache_key))
altered_manifest <- manifest
altered_manifest$design$md5[1] <- "invalid_checksum"
stopifnot(inherits(try(read_dataset(altered_manifest, altered_manifest$design$case_id[1]), silent = TRUE), "try-error"))
cat("PASS: changing q uses distinct method cache keys; checksum mismatches are rejected.\n")

registry <- load_method_registry(PROJECT_ROOT)
registry <- register_method(registry, "BrokenMethod", function(input, q, options) stop("intentional integration-test failure"))
config$methods <- "BrokenMethod"
failed <- run_experiment(manifest, config, PROJECT_ROOT, registry)
stopifnot(all(failed$metrics$status == "error"), all(is.na(failed$metrics$FDP)),
          all(failed$summary$n_error == 2L), all(failed$summary$n_success == 0L),
          all(is.na(failed$summary$FDR_mean)))
cat("PASS: method failures are diagnosed and stay NA, never converted into zero FDP.\n")
cat("Pipeline checks passed.\n")
