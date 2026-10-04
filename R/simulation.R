build_design <- function(config) {
  validate_config(config)
  grid <- expand.grid(scenario = config$scenarios, mixture = names(config$mixtures),
                      n = config$n, tau_index = seq_along(config$tau),
                      replicate = seq_len(config$repetitions), stringsAsFactors = FALSE)
  grid$tau <- config$tau[grid$tau_index]
  grid$m <- config$m
  grid$profile <- config$dgp_profile
  grid$case_id <- sprintf("%s_%s_n%d_t%02d_r%04d", grid$scenario, grid$mixture,
                          grid$n, grid$tau_index, grid$replicate)
  # Stable per-case seeds: independent of execution order and number of methods.
  grid$seed <- vapply(seq_len(nrow(grid)), function(i) {
    case <- grid[i, c("scenario", "mixture", "n", "tau", "m", "profile", "replicate")]
    h <- object_hash(list(base_seed = config$seed, case = as.list(case)))
    as.integer(strtoi(substr(h, 1, 7), base = 16) + 1L)
  }, integer(1))
  grid
}

simulate_dataset <- function(case, config) {
  n <- as.integer(case$n)
  m <- as.integer(case$m)
  with_seed(as.integer(case$seed), function() {
    X <- stats::rbinom(n, 1, config$exposure_probability)
    # A degenerate exposure is a failed dataset, never silently conditioned away.
    if (length(unique(X)) < 2L) stop("Exposure has no variation; use another seed or larger n.")
    Z <- if (case$scenario == "confounded") stats::rnorm(n) else numeric(0)
    states <- sample(c("H00", "H10", "H01", "H11"), m, replace = TRUE,
                     prob = config$mixtures[[case$mixture]])
    if (config$dgp_profile %in% c("paper2026", "paper_adjusted")) {
      sd_alpha <- sqrt(config$alpha_noise_variance / n)
      sd_beta <- sqrt(config$beta_noise_variance / n)
    } else {
      # Matches rnorm(m, mu1, kap=1), rnorm(m, mu2, psi=4) in author scripts.
      sd_alpha <- config$alpha_noise_variance
      sd_beta <- config$beta_noise_variance
    }
    alpha <- stats::rnorm(m, config$alpha_mean_scale * case$tau, sd_alpha) * (states %in% c("H10", "H11"))
    beta <- stats::rnorm(m, config$beta_mean_scale * case$tau, sd_beta) * (states %in% c("H01", "H11"))
    direct <- stats::rnorm(m, config$direct_mean, config$direct_sd)
    theta <- delta <- rep(0, m)
    if (case$scenario == "confounded") {
      theta <- stats::runif(m, 0, config$confounder_max)
      delta <- stats::runif(m, 0, config$confounder_max)
    }
    # Rows = subjects; columns = candidate mediator-outcome pathways.
    M <- outer(X, alpha) + matrix(stats::rnorm(n * m, sd = config$error_sd_m), n, m)
    if (length(Z)) M <- M + outer(Z, theta)
    eta <- sweep(M, 2, beta, "*") + outer(X, direct)
    if (length(Z)) eta <- eta + outer(Z, delta)
    Y <- if (case$scenario == "binary") {
      matrix(stats::rbinom(n * m, 1, stats::plogis(eta)), n, m)
    } else eta + matrix(stats::rnorm(n * m, sd = config$error_sd_y), n, m)
    ids <- sprintf("pathway_%04d", seq_len(m))
    colnames(M) <- colnames(Y) <- ids
    list(schema_version = 1L, meta = as.list(case),
         X = X, Z = Z, M = M, Y = Y,
         truth = data.frame(pathway_id = ids, state = states, alpha = alpha, beta = beta,
                            direct = direct, theta = theta, delta = delta,
                            is_nonnull = states == "H11", stringsAsFactors = FALSE),
         dgp = dgp_config(config))
  })
}

export_long_data <- function(data, path) {
  n <- nrow(data$M)
  m <- ncol(data$M)
  table <- data.frame(subject_id = rep(seq_len(n), m),
                      pathway_id = rep(colnames(data$M), each = n),
                      X = rep(data$X, m),
                      Z = if (length(data$Z)) rep(data$Z, m) else NA_real_,
                      M = as.vector(data$M), Y = as.vector(data$Y))
  write_csv(table, path)
  write_csv(data$truth, sub("[.]csv$", "_truth.csv", path))
  invisible(table)
}

generate_data_files <- function(config, root, force = FALSE, verbose = TRUE) {
  design <- build_design(config)
  generator_signature <- code_signature(root, file.path(root, "R", c("utils.R", "config.R", "simulation.R")))
  data_id <- paste0("data_", substr(object_hash(list(dgp = dgp_config(config), generator = generator_signature)), 1, 12))
  directory <- ensure_dir(file.path(root, "data", data_id))
  design$file <- paste0(design$case_id, ".rds")
  design$md5 <- NA_character_
  for (i in seq_len(nrow(design))) {
    path <- file.path(directory, design$file[i])
    if (force || !file.exists(path)) {
      data <- simulate_dataset(design[i, ], config)
      atomic_save_rds(data, path)
    } else {
      data <- readRDS(path)
      if (!identical(data$dgp, dgp_config(config)) || data$meta$seed != design$seed[i])
        stop("Existing data metadata mismatch: ", path)
    }
    design$md5[i] <- file_hash(path)
    if (verbose && (i == 1 || i %% 20L == 0L || i == nrow(design)))
      message(sprintf("Data saved/checked: %d/%d", i, nrow(design)))
  }
  manifest <- list(schema_version = 1L, data_id = data_id,
                   directory = directory, design = design, dgp = dgp_config(config),
                   generator_signature = generator_signature)
  atomic_save_rds(manifest, file.path(directory, "manifest.rds"))
  write_csv(design, file.path(directory, "manifest.csv"))
  if (isTRUE(config$export_first_long_csv)) {
    csv <- file.path(directory, "example_long.csv")
    if (force || !file.exists(csv)) export_long_data(readRDS(file.path(directory, design$file[1])), csv)
  }
  manifest
}

read_manifest <- function(path) {
  manifest <- readRDS(path)
  manifest$directory <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  stopifnot(manifest$schema_version == 1L, !anyDuplicated(manifest$design$case_id))
  manifest
}

read_dataset <- function(manifest, case_id, verify = TRUE) {
  idx <- match(case_id, manifest$design$case_id)
  if (is.na(idx)) stop("Unknown case_id: ", case_id)
  path <- file.path(manifest$directory, manifest$design$file[idx])
  if (verify && file_hash(path) != manifest$design$md5[idx]) stop("Data checksum mismatch: ", path)
  data <- readRDS(path)
  stopifnot(nrow(data$M) == data$meta$n, ncol(data$M) == data$meta$m,
            identical(dim(data$M), dim(data$Y)), nrow(data$truth) == data$meta$m)
  data
}
