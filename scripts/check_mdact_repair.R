PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
old_dir <- file.path(PROJECT_ROOT, "results", "run_7c638fcb10f3")
old <- readRDS(file.path(old_dir, "experiment.rds"))
manifest <- readRDS(file.path(old_dir, "provenance.rds"))$manifest
failed <- old$metrics[old$metrics$status == "error", ]
stopifnot(nrow(failed) == 1, failed$method == "MDACT")
data <- read_dataset(manifest, failed$case_id)
estimates <- estimate_coefficients(data)$estimates
fixed <- capture_conditions(function() run_mdact(make_method_input(data, estimates), failed$q))
stopifnot(!inherits(fixed$value, "method_error"))
report <- cbind(data.frame(case_id=failed$case_id, original_seed=failed$seed, previous_error=failed$error,
                          repaired_status="ok", same_saved_data=TRUE),
                as.data.frame(compute_metrics(fixed$value$reject, data$truth)))
write_csv(report, file.path(PROJECT_ROOT, "docs", "audit", "mdact_repair.csv"))
print(report, row.names=FALSE)
cat("PASS: original failed dataset was rerun without changing its data or seed.\n")
