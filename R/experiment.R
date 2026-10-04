compute_metrics <- function(reject, truth) {
  stopifnot(is.logical(reject), length(reject) == nrow(truth), !anyNA(reject))
  true <- truth$is_nonnull
  R <- sum(reject)
  V <- sum(reject & !true)
  S <- sum(reject & true)
  alternatives <- sum(true)
  list(discoveries = R, false_discoveries = V, true_discoveries = S,
       alternatives = alternatives, FDP = V / max(R, 1L),
       power = if (alternatives > 0L) S / alternatives else NA_real_)
}

method_fingerprint <- function(entry, root) {
  files <- file.path(root, c("R/utils.R", "R/methods/builtin.R", "R/methods/mlfdr_em.R", "R/methods/registry.R", "vendor/mdact_core.R"))
  list(name = entry$name, version = entry$version, code = code_signature(root, files),
       body = deparse(body(entry$run)), formals = formals(entry$run), source_hash = entry$source_hash,
       packages = package_versions(unique(c(entry$packages, "HDMT", "MLFDR"))))
}

analyse_case <- function(index, manifest, config, root, registry) {
  case <- manifest$design[index, ]
  data <- read_dataset(manifest, case$case_id)
  estimate_key <- object_hash(list(data_md5 = case$md5,
                                  code = code_signature(root, file.path(root, "R", c("utils.R", "estimation.R"))),
                                  R = R.version.string))
  estimates_path <- file.path(root, "data", "estimates", paste0(estimate_key, ".rds"))
  estimation <- capture_conditions(function() {
    if (file.exists(estimates_path)) readRDS(estimates_path) else {
      fit <- estimate_coefficients(data)
      atomic_save_rds(fit, estimates_path)
      fit
    }
  })
  rows <- list()
  for (method in config$methods) {
    entry <- registry[[method]]
    opts <- config$method_options[[method]]
    if (is.null(opts)) opts <- list()
    method_seed <- as.integer(strtoi(substr(object_hash(list(case$seed, method)), 1, 7), 16) + 1L)
    cache_key <- object_hash(list(estimate_key = estimate_key, q = config$q, options = opts,
                                 method_seed = method_seed, fingerprint = method_fingerprint(entry, root)))
    cache_path <- file.path(root, "results", "method_cache", method, paste0(cache_key, ".rds"))
    cached <- file.exists(cache_path)
    if (cached) {
      outcome <- readRDS(cache_path)
    } else if (inherits(estimation$value, "method_error")) {
      outcome <- list(status = "error", error = paste("Estimation:", estimation$value$message),
                      warnings = estimation$warnings, elapsed = estimation$elapsed, result = NULL)
    } else {
      input <- make_method_input(data, estimation$value$estimates)
      captured <- capture_conditions(function() with_seed(method_seed, function() {
        validate_method_result(entry$run(input, config$q, opts), ncol(data$M))
      }))
      failed <- inherits(captured$value, "method_error")
      outcome <- list(status = if (failed) "error" else "ok",
                      error = if (failed) captured$value$message else "",
                      warnings = unique(c(estimation$value$warnings, captured$warnings)),
                      elapsed = captured$elapsed, result = if (failed) NULL else captured$value)
      # Failed method calls are not cached; fixes/retries remain possible.
      if (!failed) atomic_save_rds(outcome, cache_path)
    }
    metrics <- if (outcome$status == "ok") compute_metrics(outcome$result$reject, data$truth) else {
      list(discoveries = NA_integer_, false_discoveries = NA_integer_, true_discoveries = NA_integer_,
           alternatives = sum(data$truth$is_nonnull), FDP = NA_real_, power = NA_real_)
    }
    rows[[method]] <- cbind(case[, c("case_id", "scenario", "mixture", "profile", "n", "m", "tau", "replicate", "seed")],
                             data.frame(method = method, q = config$q, status = outcome$status,
                                        cached = cached, elapsed_seconds = outcome$elapsed,
                                        warning = paste(outcome$warnings, collapse = " | "), error = outcome$error,
                                        cache_key = cache_key, stringsAsFactors = FALSE),
                             as.data.frame(metrics))
  }
  do.call(rbind, rows)
}

run_experiment <- function(manifest, config, root, registry = load_method_registry(root)) {
  validate_config(config)
  if (!identical(manifest$dgp, dgp_config(config))) stop("Config and saved data differ. Generate/read the matching manifest.")
  missing <- setdiff(config$methods, names(registry))
  if (length(missing)) stop("Unregistered methods: ", paste(missing, collapse = ", "))
  required <- unique(unlist(lapply(registry[config$methods], function(e) e$packages)))
  versions <- package_versions(required)
  if (anyNA(versions$version)) stop("Missing method dependencies: ", paste(versions$package[is.na(versions$version)], collapse = ", "))
  fingerprint <- lapply(registry[config$methods], method_fingerprint, root = root)
  run_id <- paste0("run_", substr(object_hash(list(data_id = manifest$data_id, config = config[c("q", "methods", "method_options")],
                                                methods = fingerprint,
                                                engine = code_signature(root, file.path(root, "R", c("experiment.R", "reporting.R"))))), 1, 12))
  output_dir <- ensure_dir(file.path(root, "results", run_id))
  atomic_save_rds(config, file.path(output_dir, "config.rds"))
  dput(config, file = file.path(output_dir, "config.R"))
  write_csv(registry_info(registry[config$methods]), file.path(output_dir, "methods.csv"))
  write_csv(package_versions(unique(c(required, "ggplot2", "IRkernel"))), file.path(output_dir, "package_versions.csv"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"), useBytes = TRUE)
  atomic_save_rds(list(manifest = manifest, fingerprint = fingerprint, code = code_signature(root)),
                  file.path(output_dir, "provenance.rds"))
  indices <- seq_len(nrow(manifest$design))
  worker <- function(index) analyse_case(index, manifest, config, root, registry)
  # PSOCK works on Windows; mclapply(mc.cores > 1) does not.
  workers <- min(config$workers, length(indices))
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(cluster, "root", envir = environment())
    parallel::clusterEvalQ(cluster, {
      PROJECT_ROOT <- root
      source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
      NULL
    })
  }
  rows <- list()
  chunks <- split(indices, ceiling(seq_along(indices) / 12L))
  for (chunk in chunks) {
    values <- if (workers > 1L) parallel::parLapplyLB(cluster, chunk, worker) else lapply(chunk, worker)
    rows <- c(rows, values)
    current <- do.call(rbind, rows)
    write_csv(current, file.path(output_dir, "replicate_metrics.csv"))
    message(sprintf("Analysed %d/%d saved datasets; method errors: %d",
                    length(rows), length(indices), sum(current$status == "error")))
  }
  metrics <- do.call(rbind, rows)
  rownames(metrics) <- NULL
  summary <- summarise_metrics(metrics)
  write_csv(summary, file.path(output_dir, "summary.csv"))
  diagnostics <- metrics[, c("case_id", "method", "status", "warning", "error", "elapsed_seconds", "cached")]
  write_csv(diagnostics, file.path(output_dir, "diagnostics.csv"))
  atomic_save_rds(list(metrics = metrics, summary = summary, data_id = manifest$data_id, run_id = run_id),
                  file.path(output_dir, "experiment.rds"))
  list(metrics = metrics, summary = summary, diagnostics = diagnostics,
       output_dir = output_dir, run_id = run_id, config = config, manifest = manifest)
}

get_case_results <- function(experiment, case_id, root) {
  rows <- experiment$metrics[experiment$metrics$case_id == case_id, , drop = FALSE]
  results <- lapply(seq_len(nrow(rows)), function(i) {
    path <- file.path(root, "results", "method_cache", rows$method[i], paste0(rows$cache_key[i], ".rds"))
    if (rows$status[i] != "ok" || !file.exists(path)) return(NULL)
    readRDS(path)$result
  })
  setNames(results, rows$method)
}
