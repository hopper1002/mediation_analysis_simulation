# From the project root: Rscript --vanilla scripts/run_experiment.R demo
args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args)) args[[1]] else "demo"
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash = "/")
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
check_dependencies()
config <- default_config(preset)
print(estimate_storage(config))
manifest <- generate_data_files(config, PROJECT_ROOT)
manifest <- read_manifest(file.path(manifest$directory, "manifest.rds"))
experiment <- run_experiment(manifest, config, PROJECT_ROOT)
save_figures(experiment)
if (config$preset == "paper_comparison") {
  save_comparison_artifacts(experiment, PROJECT_ROOT)
  atomic_save_rds(list(output_dir = experiment$output_dir,
                       manifest_path = file.path(manifest$directory, "manifest.rds"), config = config),
                  file.path(PROJECT_ROOT, "docs", "current_comparison.rds"))
}
cat("Results:", experiment$output_dir, "\n")
