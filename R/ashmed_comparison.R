# Comparison orchestration only: consumes saved manifests; never generates data.

ashmed_saved_sources <- function() {
  list(paper_adjusted = list(label = "Dataset 1: paper-adjusted", data_id = "data_f87521c2b4eb",
                            baseline_run = "run_2836e6c08e18"),
       teaching_fixed = list(label = "Dataset 2: teaching-fixed", data_id = "data_c054ed713a76",
                            baseline_run = "run_af52854e0ae0"))
}

read_ashmed_inputs <- function(root, sources = ashmed_saved_sources(),
                               method_options = list(ASHMED = list()), workers = 4L) {
  stopifnot("ASHMED" %in% names(method_options), !anyDuplicated(names(method_options)))
  if (any(names(method_options) %in% c("HDMT", "MDACT", "MLFDR")))
    stop("Keep the saved baseline method options unchanged; configure additional methods only.")
  inputs <- lapply(sources, function(spec) {
    manifest_path <- file.path(root, "data", spec$data_id, "manifest.rds")
    baseline_path <- file.path(root, "results", spec$baseline_run)
    required <- c(manifest_path, file.path(baseline_path, c("config.rds", "replicate_metrics.csv")))
    if (!all(file.exists(required)))
      stop("Missing saved input files. This notebook only reads data; update the paths to existing manifests: ",
           paste(required[!file.exists(required)], collapse = ", "))
    manifest <- read_manifest(manifest_path)
    config <- readRDS(file.path(baseline_path, "config.rds"))
    if (!identical(manifest$dgp, dgp_config(config)))
      stop("Saved config does not match the saved manifest.")
    if (!identical(config$methods, c("HDMT", "MDACT", "MLFDR")))
      stop("Expected the saved three-method baseline configuration.")
    baseline <- utils::read.csv(file.path(baseline_path, "replicate_metrics.csv"), stringsAsFactors = FALSE)
    if (any(baseline$status != "ok") || nrow(baseline) != 3L * nrow(manifest$design))
      stop("Baseline is incomplete or contains failed methods.")
    cache_paths <- file.path(root, "results/method_cache", baseline$method, paste0(baseline$cache_key, ".rds"))
    if (!all(file.exists(cache_paths))) stop("Existing baseline method caches are missing; restore them first.")
    config$methods <- c(config$methods, names(method_options))
    for (method in names(method_options)) config$method_options[[method]] <- method_options[[method]]
    config$workers <- as.integer(workers)
    validate_config(config)
    raw_paths <- file.path(manifest$directory, manifest$design$file)
    if (!all(file.exists(raw_paths))) stop("Some saved raw datasets are missing.")
    list(spec = spec, config = config, manifest = manifest, manifest_path = manifest_path,
         manifest_md5 = file_hash(manifest_path), raw_paths = raw_paths,
         original_mtime = file.info(raw_paths)$mtime,
         baseline = baseline, baseline_path = baseline_path,
         baseline_md5 = file_hash(file.path(baseline_path, "replicate_metrics.csv")))
  })
  inputs
}

audit_ashmed_baseline <- function(experiment, input) {
  columns <- c("discoveries", "false_discoveries", "true_discoveries", "alternatives", "FDP", "power")
  baseline <- input$baseline
  new <- experiment$metrics[experiment$metrics$method %in% unique(baseline$method), ]
  key <- function(x) paste(x$case_id, x$method, sep = "::")
  stopifnot(!anyDuplicated(key(new)), setequal(key(new), key(baseline)))
  new <- new[match(key(baseline), key(new)), ]
  result <- do.call(rbind, lapply(columns, function(name) {
    delta <- abs(new[[name]] - baseline[[name]])
    data.frame(metric = name, n_compared = length(delta), max_absolute_difference = max(delta),
               unchanged = all(delta <= 1e-12))
  }))
  if (!all(result$unchanged) || any(new$status != "ok") || !all(new$cached))
    stop("Original methods changed or were recomputed. Inspect the baseline audit before reporting.")
  if (file_hash(file.path(input$baseline_path, "replicate_metrics.csv")) != input$baseline_md5)
    stop("Original baseline output was modified.")
  result
}

audit_ashmed_data <- function(input) {
  md5 <- unname(tools::md5sum(input$raw_paths))
  result <- data.frame(data_id = input$manifest$data_id, datasets = length(md5),
                       checksum_matches = sum(md5 == input$manifest$design$md5),
                       mtime_unchanged = all(file.info(input$raw_paths)$mtime == input$original_mtime),
                       manifest_unchanged = file_hash(input$manifest_path) == input$manifest_md5)
  stopifnot(result$checksum_matches == result$datasets, result$mtime_unchanged, result$manifest_unchanged)
  result
}

collect_ashmed_fit_diagnostics <- function(experiment, root) {
  rows <- experiment$metrics[experiment$metrics$method == "ASHMED" & experiment$metrics$status == "ok", ]
  output <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    result <- readRDS(file.path(root, "results/method_cache/ASHMED", paste0(row$cache_key, ".rds")))$result
    fit <- result$diagnostics
    data.frame(case_id = row$case_id, scenario = row$scenario, mixture = row$mixture, tau = row$tau,
               converged = fit$converged, iterations = fit$iterations, dual_gap = fit$dual_gap,
               loglik = fit$loglik, pi00 = fit$pi[1], pi10 = fit$pi[2], pi01 = fit$pi[3], pi11 = fit$pi[4],
               discoveries = sum(result$reject), selected_mean_lfdr = fit$selected_mean_lfdr,
               n_fit = fit$n_fit, n_components = nrow(result$components), likelihood = fit$likelihood)
  })
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}

plot_ashmed_comparison <- function(summary, scenario) {
  colors <- c(HDMT = "#2864A0", MDACT = "#13856B", MLFDR = "#C1486A", ASHMED = "#9C6900")
  extra <- setdiff(unique(summary$method), names(colors))
  if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  # Public ggplot API replaces only the palette; all original panels/metrics stay the same.
  suppressMessages(plot_paper_comparison(summary, scenario) +
                     ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors))
}

export_ashmed_example <- function(experiment, root) {
  design <- experiment$manifest$design
  idx <- which(design$scenario == "linear" & design$mixture == "dense" &
                 design$tau == 1.5 & design$replicate == 1L)[1]
  if (is.na(idx)) idx <- which(design$scenario == "linear")[1]
  if (is.na(idx)) idx <- 1L
  case_id <- design$case_id[idx]
  data <- read_dataset(experiment$manifest, case_id)
  estimates <- estimate_coefficients(data)$estimates
  result <- get_case_results(experiment, case_id, root)$ASHMED
  # Truth is joined AFTER fitting and selection, only for illustration/evaluation.
  joined <- cbind(estimates[, c("pathway_id", "alpha_hat", "beta_hat", "se_alpha", "se_beta")],
                   result$posterior[, -1],
                   truth_state = data$truth$state, truth_alpha = data$truth$alpha, truth_beta = data$truth$beta)
  write_csv(joined, file.path(experiment$output_dir, "ashmed_example_posterior.csv"))
  write_csv(result$components, file.path(experiment$output_dir, "ashmed_example_components.csv"))
  write_csv(result$diagnostics$trace, file.path(experiment$output_dir, "ashmed_example_optimizer.csv"))
  list(case_id = case_id, posterior = joined, result = result)
}

run_ashmed_suite <- function(inputs, root, registry = load_method_registry(root)) {
  output <- list()
  for (label in names(inputs)) {
    input <- inputs[[label]]
    message("Four-method comparison on saved data: ", label)
    experiment <- run_experiment(input$manifest, input$config, root, registry)
    baseline_audit <- audit_ashmed_baseline(experiment, input)
    data_audit <- audit_ashmed_data(input)
    write_csv(baseline_audit, file.path(experiment$output_dir, "baseline_audit.csv"))
    write_csv(data_audit, file.path(experiment$output_dir, "data_audit.csv"))
    paired <- summarise_paired_metrics(experiment$metrics, left = "ASHMED",
                                       right = setdiff(experiment$config$methods, "ASHMED"))
    fits <- collect_ashmed_fit_diagnostics(experiment, root)
    write_csv(paired, file.path(experiment$output_dir, "paired_differences.csv"))
    write_csv(fits, file.path(experiment$output_dir, "ashmed_fit_diagnostics.csv"))
    # Preserve the established MLFDR diagnostics alongside the new diagnostics.
    write_csv(collect_fit_diagnostics(experiment, root), file.path(experiment$output_dir, "fit_diagnostics.csv"))
    directory <- ensure_dir(file.path(experiment$output_dir, "figures"))
    for (scenario in unique(experiment$summary$scenario))
      ggplot2::ggsave(file.path(directory, paste0(scenario, "_comparison.png")),
                      plot_ashmed_comparison(experiment$summary, scenario), width = 10, height = 7, dpi = 160, bg = "white")
    example <- export_ashmed_example(experiment, root)
    output[[label]] <- list(experiment = experiment, paired = paired, fits = fits,
                            example = example, baseline_audit = baseline_audit, data_audit = data_audit)
  }
  pointer <- lapply(output, function(x) list(run_dir = x$experiment$output_dir,
                                           manifest_path = file.path(x$experiment$manifest$directory, "manifest.rds")))
  atomic_save_rds(pointer, file.path(root, "docs/current_ashmed.rds"))
  writeLines(unlist(lapply(names(pointer), function(name) c(paste0("dataset=", name),
               paste0("run_dir=", pointer[[name]]$run_dir), paste0("manifest_path=", pointer[[name]]$manifest_path)))),
             file.path(root, "docs/current_ashmed.txt"), useBytes = TRUE)
  output
}
