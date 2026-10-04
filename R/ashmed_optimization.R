# Candidate development and held-out evaluation on already saved repeats.
optimization_directory <- function(root) ensure_dir(file.path(root, "results/ashmed_optimization_20261005"))

optimization_inputs <- function(root) {
  specs <- list(paper_adjusted = c("data_f87521c2b4eb", "run_76428dee9701"),
                teaching_fixed = c("data_c054ed713a76", "run_17d85e6975d0"))
  lapply(specs, function(spec) {
    manifest_path <- file.path(root, "data", spec[1], "manifest.rds")
    original_dir <- file.path(root, "results", spec[2])
    manifest <- read_manifest(manifest_path)
    config <- readRDS(file.path(original_dir, "config.rds"))
    reference <- readRDS(file.path(original_dir, "experiment.rds"))$metrics
    stopifnot(identical(manifest$dgp, dgp_config(config)), all(reference$status == "ok"))
    list(manifest = manifest, config = config, reference = reference, original_dir = original_dir,
         manifest_path = manifest_path, manifest_md5 = file_hash(manifest_path))
  })
}

run_adaptive_trial <- function(inputs, root, options, stage, replicates, arm, workers = 4L) {
  stopifnot(length(replicates) > 1L, !anyDuplicated(replicates), all(replicates %in% 1:40))
  results <- list()
  registry <- load_method_registry(root)
  canonical <- environment(registry$ASHMED_Adaptive$run)$adaptive_options(options)
  canonical$base <- NULL
  for (dataset in names(inputs)) {
    input <- inputs[[dataset]]
    manifest <- input$manifest
    manifest$source_data_id <- manifest$data_id
    manifest$design <- manifest$design[manifest$design$replicate %in% replicates, , drop = FALSE]
    # Stage in data_id prevents the old engine's run hash from confusing
    # distinct subsets of the same manifest. Raw directory/MD5 stay unchanged.
    manifest$data_id <- paste0(manifest$data_id, "_", stage, "_r", paste(replicates, collapse = "-"))
    config <- input$config
    config$methods <- "ASHMED_Adaptive"
    config$method_options <- list(ASHMED_Adaptive = canonical)
    config$workers <- as.integer(workers)
    config$comparison_stage <- stage
    config$analysed_replicates <- replicates
    experiment <- run_experiment(manifest, config, root, registry)
    references <- input$reference[input$reference$case_id %in% manifest$design$case_id, ]
    references$cached <- TRUE
    combined <- rbind(references, experiment$metrics)
    combined$dataset <- dataset
    combined$stage <- stage
    combined$arm <- arm
    summary <- summarise_metrics(combined)
    summary$dataset <- dataset; summary$stage <- stage; summary$arm <- arm
    results[[dataset]] <- list(experiment = experiment, combined = combined, summary = summary)
  }
  dir <- optimization_directory(root)
  write_csv(do.call(rbind, lapply(results, `[[`, "combined")), file.path(dir, paste0(stage, "_", arm, "_metrics.csv")))
  write_csv(do.call(rbind, lapply(results, `[[`, "summary")), file.path(dir, paste0(stage, "_", arm, "_summary.csv")))
  atomic_save_rds(list(results = results, options = options, stage = stage, replicates = replicates, arm = arm),
                  file.path(dir, paste0(stage, "_", arm, ".rds")))
  results
}

score_adaptive_trial <- function(results, arm) {
  do.call(rbind, lapply(names(results), function(dataset) {
    s <- results[[dataset]]$summary
    a <- s[s$method == "ASHMED_Adaptive", ]
    original <- s[s$method == "ASHMED", ]
    data.frame(arm = arm, dataset = dataset, mean_FDR = mean(a$FDR_mean),
      max_FDR = max(a$FDR_mean), mean_power = mean(a$power_mean),
      original_power = mean(original$power_mean), gain = mean(a$power_mean) - mean(original$power_mean),
      errors = sum(a$n_error), warnings = sum(a$n_warning),
      passes = mean(a$FDR_mean) <= .06 && max(a$FDR_mean) <= .15 &&
        mean(a$power_mean) > mean(original$power_mean) && sum(a$n_error) == 0)
  }))
}

plot_adaptive_comparison <- function(summary, scenario) {
  colors <- c(HDMT = "#2864A0", MDACT = "#13856B", MLFDR = "#C1486A",
              ASHMED = "#9C6900", ASHMED_Adaptive = "#6D45A0")
  suppressMessages(plot_paper_comparison(summary, scenario) +
    ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors))
}

adaptive_macro_summary <- function(metrics) {
  do.call(rbind, lapply(unique(metrics$method), function(method) {
    d <- metrics[metrics$method == method, ]
    stopifnot(all(d$status == "ok"))
    row <- data.frame(method = method, conditions = 18L, repetitions = length(unique(d$replicate)))
    for (metric in c("FDP", "power")) {
      # Equally weighted fixed conditions; uncertainty across replicate blocks.
      blocks <- tapply(d[[metric]], d$replicate, mean)
      se <- sd(blocks) / sqrt(length(blocks))
      half <- qt(.975, length(blocks) - 1) * se
      prefix <- if (metric == "FDP") "FDR" else "power"
      row[[paste0(prefix, "_mean")]] <- mean(blocks)
      row[[paste0(prefix, "_mcse")]] <- se
      row[[paste0(prefix, "_lo")]] <- max(0, mean(blocks) - half)
      row[[paste0(prefix, "_hi")]] <- min(1, mean(blocks) + half)
    }
    row
  }))
}

collect_adaptive_diagnostics <- function(experiment, root) {
  do.call(rbind, lapply(seq_len(nrow(experiment$metrics)), function(i) {
    row <- experiment$metrics[i, ]
    if (row$status != "ok") stop("Adaptive method failed: ", row$case_id)
    result <- readRDS(file.path(root, "results/method_cache/ASHMED_Adaptive",
                               paste0(row$cache_key, ".rds")))$result
    d <- result$diagnostics
    data.frame(case_id = row$case_id, scenario = row$scenario, mixture = row$mixture, tau = row$tau,
      converged = d$converged, iterations = d$iterations, dual_gap = d$dual_gap,
      location_converged = all(vapply(d$locations, function(x) x$converged, logical(1))),
      mean_alpha = d$locations$alpha$mean, mean_beta = d$locations$beta$mean,
      fallback = d$fallback, n_components = nrow(result$components), loglik = d$loglik,
      pi00 = d$pi[1], pi10 = d$pi[2], pi01 = d$pi[3], pi11 = d$pi[4],
      selected_mean_lfdr = d$selected_mean_lfdr, discoveries = sum(result$reject))
  }))
}

audit_protected_adaptive_files <- function(root) {
  before <- readRDS(file.path(optimization_directory(root), "protected_files.rds"))
  after <- unname(tools::md5sum(file.path(root, names(before))))
  result <- data.frame(file = names(before), unchanged = unname(before) == after)
  stopifnot(all(result$unchanged))
  result
}

run_final_adaptive_comparison <- function(root, workers = 4L) {
  campaign <- optimization_directory(root)
  lock_path <- file.path(campaign, "locked_method.rds")
  if (!file.exists(lock_path)) stop("Finish development and confirmation before final evaluation.")
  lock <- readRDS(lock_path)
  selection <- lock$selection
  registry <- load_method_registry(root)
  stopifnot(identical(method_fingerprint(registry$ASHMED_Adaptive, root), selection$fingerprint),
            all(lock$confirmation$passes), identical(selection$validation_replicates, 13:40))
  inputs <- optimization_inputs(root)
  trial <- run_adaptive_trial(inputs, root, selection$options, "validation",
                              selection$validation_replicates, selection$arm, workers)
  output <- list()
  for (label in names(trial)) {
    x <- trial[[label]]
    dir <- ensure_dir(file.path(campaign, paste0("validation_", label)))
    write_csv(x$combined, file.path(dir, "replicate_metrics.csv"))
    write_csv(x$summary, file.path(dir, "summary.csv"))
    macro <- adaptive_macro_summary(x$combined)
    write_csv(macro, file.path(dir, "macro_summary.csv"))
    main <- x$summary[x$summary$method != "ASHMED", ]
    write_csv(main, file.path(dir, "four_method_summary.csv"))
    paired <- summarise_paired_metrics(x$combined, left = "ASHMED_Adaptive",
                                       right = c("HDMT", "MDACT", "MLFDR", "ASHMED"))
    write_csv(paired, file.path(dir, "paired_differences.csv"))
    fits <- collect_adaptive_diagnostics(x$experiment, root)
    write_csv(fits, file.path(dir, "adaptive_fit_diagnostics.csv"))
    stopifnot(all(fits$converged), all(fits$location_converged))
    references <- inputs[[label]]$reference
    references <- references[references$case_id %in% x$experiment$manifest$design$case_id, ]
    actual <- x$combined[x$combined$method != "ASHMED_Adaptive", ]
    key <- function(d) paste(d$case_id, d$method, sep = "::")
    actual <- actual[match(key(references), key(actual)), ]
    audit <- do.call(rbind, lapply(c("discoveries", "false_discoveries", "true_discoveries",
                                   "alternatives", "FDP", "power"), function(k) {
      delta <- abs(actual[[k]] - references[[k]])
      data.frame(metric = k, n_compared = length(delta), max_difference = max(delta),
                 unchanged = all(delta < 1e-12))
    }))
    stopifnot(all(audit$unchanged), all(actual$cached))
    write_csv(audit, file.path(dir, "reference_audit.csv"))
    manifest <- inputs[[label]]$manifest
    hashes <- unname(tools::md5sum(file.path(manifest$directory, manifest$design$file)))
    data_audit <- data.frame(dataset = label, saved_datasets = length(hashes),
      checksum_matches = sum(hashes == manifest$design$md5),
      manifest_unchanged = file_hash(inputs[[label]]$manifest_path) == inputs[[label]]$manifest_md5)
    stopifnot(data_audit$checksum_matches == 720L, data_audit$manifest_unchanged)
    write_csv(data_audit, file.path(dir, "data_audit.csv"))
    figures <- ensure_dir(file.path(dir, "figures"))
    for (scenario in c("linear", "confounded", "binary"))
      ggplot2::ggsave(file.path(figures, paste0(scenario, "_comparison.png")),
        plot_adaptive_comparison(main, scenario), width = 10, height = 7, dpi = 160, bg = "white")
    original <- x$summary[x$summary$method %in% c("ASHMED", "ASHMED_Adaptive"), ]
    ggplot2::ggsave(file.path(figures, "original_vs_adaptive.png"),
      plot_adaptive_comparison(original, "linear"), width = 10, height = 7, dpi = 160, bg = "white")
    design <- x$experiment$manifest$design
    id <- design$case_id[design$scenario == "linear" & design$mixture == "dense" &
                          design$tau == 1.5 & design$replicate == 13L][1]
    data <- read_dataset(x$experiment$manifest, id)
    estimates <- estimate_coefficients(data)$estimates
    result <- get_case_results(x$experiment, id, root)$ASHMED_Adaptive
    # Truth is attached only after the empirical Bayes fit and decision.
    posterior <- cbind(estimates[, c("pathway_id", "alpha_hat", "beta_hat", "se_alpha", "se_beta")],
      result$posterior[, -1], truth_state = data$truth$state,
      truth_alpha = data$truth$alpha, truth_beta = data$truth$beta)
    write_csv(posterior, file.path(dir, "adaptive_example_posterior.csv"))
    write_csv(result$components, file.path(dir, "adaptive_example_components.csv"))
    write_csv(result$diagnostics$trace, file.path(dir, "adaptive_example_optimizer.csv"))
    output[[label]] <- c(x, list(directory = dir, main = main, macro = macro, paired = paired, fits = fits,
      reference_audit = audit, data_audit = data_audit,
      example = list(case_id = id, posterior = posterior, result = result)))
    atomic_save_rds(output[[label]], file.path(dir, "comparison.rds"))
  }
  protected <- audit_protected_adaptive_files(root)
  write_csv(protected, file.path(campaign, "protected_file_audit.csv"))
  pointer <- lapply(output, function(x) list(directory = x$directory,
    adaptive_run_dir = x$experiment$output_dir,
    manifest_path = file.path(x$experiment$manifest$directory, "manifest.rds")))
  atomic_save_rds(pointer, file.path(root, "docs/current_ashmed_adaptive.rds"))
  writeLines(unlist(lapply(names(pointer), function(n) c(paste0("dataset=", n),
    paste0("directory=", pointer[[n]]$directory), paste0("adaptive_run_dir=", pointer[[n]]$adaptive_run_dir)))),
    file.path(root, "docs/current_ashmed_adaptive.txt"), useBytes = TRUE)
  output
}
