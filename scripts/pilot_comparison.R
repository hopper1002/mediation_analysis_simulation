# Development pilot; independent base seed from the delivered comparison.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
cfg <- default_config("paper_comparison")
cfg$seed <- 20261005L
cfg$dgp_profile <- "paper2026" # Historical pilot tag; alpha_mean_scale=.35 is still adjusted.
cfg$repetitions <- 4L
cfg$n <- c(100L, 300L)
cfg$tau <- c(0.7, 1.3, 1.9)
cfg$scenarios <- c("linear", "binary")
manifest <- generate_data_files(cfg, PROJECT_ROOT)
result <- run_experiment(manifest, cfg, PROJECT_ROOT)
save_figures(result)
cat("PILOT OUTPUT:", result$output_dir, "\n")
print(result$summary[, c("scenario", "mixture", "n", "tau", "method", "FDR_mean", "power_mean", "n_error", "n_warning")], row.names = FALSE)
