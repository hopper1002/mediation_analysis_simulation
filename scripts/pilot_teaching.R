PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
config <- default_config("teaching")
config$seed <- 20261008L
config$repetitions <- 3L
manifest <- generate_data_files(config, PROJECT_ROOT)
experiment <- run_experiment(manifest, config, PROJECT_ROOT)
artifacts <- save_comparison_artifacts(experiment, PROJECT_ROOT)
atomic_save_rds(list(output_dir=experiment$output_dir,
                     manifest_path=file.path(manifest$directory, "manifest.rds"), config=config),
                file.path(PROJECT_ROOT, "docs", "teaching_pilot.rds"))
print(experiment$summary[,c("scenario","mixture","tau","method","n_error","FDR_mean","power_mean")], row.names=FALSE)
cat("Errors:", sum(experiment$metrics$status != "ok"), "\n")
cat("EM converged:", sum(artifacts$fits$converged), "/", nrow(artifacts$fits), "\n")
cat("Run:", experiment$output_dir, "\n")
