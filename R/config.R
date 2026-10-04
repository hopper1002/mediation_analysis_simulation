default_config <- function(preset = c("demo", "paper_grid", "source_demo", "paper_comparison", "teaching")) {
  preset <- match.arg(preset)
  config <- list(
    schema_version = 1L,
    preset = preset,
    seed = 20261004L,
    m = 1000L,
    n = c(100L, 300L),
    tau = c(0.1, 1.0, 1.9),
    repetitions = 5L,
    scenarios = c("linear", "confounded", "binary"),
    mixtures = list(sparse = c(H00 = 0.88, H10 = 0.05, H01 = 0.05, H11 = 0.02),
                    dense = c(H00 = 0.4, H10 = 0.2, H01 = 0.2, H11 = 0.2)),
    dgp_profile = "paper2026",
    exposure_probability = 0.1,
    alpha_mean_scale = 0.05,
    beta_mean_scale = -0.5,
    alpha_noise_variance = 1,
    beta_noise_variance = 4,
    direct_mean = 1,
    direct_sd = sqrt(0.5),
    error_sd_m = 1,
    error_sd_y = 1,
    confounder_max = 0.5,
    q = 0.05,
    methods = c("HDMT", "MDACT", "MLFDR"),
    method_options = list(HDMT = list(exact = 0L), MDACT = list(),
                          MLFDR = list(eps = 0.01, twostep = FALSE, verbose = FALSE)),
    workers = 2L,
    export_first_long_csv = TRUE
  )
  if (preset == "paper_grid") {
    config$tau <- seq(0.1, 1.9, by = 0.2)
    config$repetitions <- 250L
  }
  if (preset == "source_demo") {
    config$dgp_profile <- "source_code"
    config$direct_mean <- 0.5
    config$direct_sd <- 1
  }
  if (preset %in% c("paper_comparison", "teaching")) {
    config$seed <- 20261007L
    config$dgp_profile <- "paper_adjusted"
    config$alpha_mean_scale <- 0.35
    config$n <- 300L
    config$tau <- c(1.3, 1.5, 1.9)
    config$repetitions <- 40L
    config$workers <- 4L
    config$method_options$MLFDR <- list(engine = "paper_em", eps = 1e-4, max_iter = 2000L, n_starts = 3L)
  }
  if (preset == "teaching") {
    config$seed <- 20261010L
    config$dgp_profile <- "teaching_fixed"
    config$exposure_probability <- 0.5
    config$alpha_mean_scale <- 0.2
    config$beta_mean_scale <- 0.5
    config$alpha_noise_variance <- config$beta_noise_variance <- 0
    config$direct_mean <- 0.2
    config$direct_sd <- 0
    config$confounder_max <- 0.3
  }
  config
}

validate_config <- function(config) {
  is_int <- function(x) is.numeric(x) && length(x) > 0 && all(is.finite(x) & x == floor(x))
  stopifnot(is_int(config$m), length(config$m) == 1, config$m >= 10,
            is_int(config$n), all(config$n >= 10), !anyDuplicated(config$n),
            is_int(config$repetitions), length(config$repetitions) == 1, config$repetitions >= 1,
            is_int(config$seed), length(config$seed) == 1, config$seed >= 0,
            is_int(config$workers), length(config$workers) == 1, config$workers >= 1,
            is.numeric(config$tau), length(config$tau) > 0, all(is.finite(config$tau)),
            !anyDuplicated(config$tau), length(config$q) == 1, config$q > 0, config$q < 1,
            all(config$scenarios %in% c("linear", "confounded", "binary")),
            length(config$scenarios) > 0, !anyDuplicated(config$scenarios),
            config$dgp_profile %in% c("paper2026", "paper_adjusted", "source_code", "teaching_fixed"),
            length(config$methods) > 0, !anyDuplicated(config$methods),
            all(grepl("^[A-Za-z][A-Za-z0-9_]*$", config$methods)),
            length(config$exposure_probability) == 1,
            config$exposure_probability > 0, config$exposure_probability < 1)
  stopifnot(length(config$mixtures) > 0, !is.null(names(config$mixtures)),
            all(grepl("^[A-Za-z][A-Za-z0-9_]*$", names(config$mixtures))),
            !anyDuplicated(names(config$mixtures)))
  for (pi in config$mixtures) {
    stopifnot(identical(names(pi), c("H00", "H10", "H01", "H11")),
              all(is.finite(pi)), all(pi >= 0), abs(sum(pi) - 1) < 1e-10)
  }
  for (key in c("alpha_noise_variance", "beta_noise_variance", "direct_sd")) {
    stopifnot(length(config[[key]]) == 1, is.finite(config[[key]]))
    if (config$dgp_profile == "teaching_fixed") stopifnot(config[[key]] == 0)
    else stopifnot(config[[key]] > 0)
  }
  for (key in c("error_sd_m", "error_sd_y"))
    stopifnot(length(config[[key]]) == 1, is.finite(config[[key]]), config[[key]] > 0)
  if (config$dgp_profile == "teaching_fixed")
    stopifnot(all(config$tau > 0), config$alpha_mean_scale != 0, config$beta_mean_scale != 0)
  for (key in c("alpha_mean_scale", "beta_mean_scale", "direct_mean", "confounder_max"))
    stopifnot(length(config[[key]]) == 1, is.finite(config[[key]]))
  stopifnot(config$confounder_max >= 0)
  invisible(config)
}

# Data settings exclude inference settings: changing q or adding methods preserves data.
dgp_config <- function(config) {
  keys <- c("schema_version", "seed", "m", "n", "tau", "repetitions", "scenarios", "mixtures",
            "dgp_profile", "exposure_probability", "alpha_mean_scale", "beta_mean_scale",
            "alpha_noise_variance", "beta_noise_variance", "direct_mean", "direct_sd",
            "error_sd_m", "error_sd_y", "confounder_max")
  config[keys]
}

design_overview <- function(config) {
  grid <- expand.grid(scenario = config$scenarios, mixture = names(config$mixtures),
                      n = config$n, tau = config$tau, stringsAsFactors = FALSE)
  grid$m <- config$m
  grid$repetitions <- config$repetitions
  grid$profile <- config$dgp_profile
  grid
}

estimate_storage <- function(config) {
  design <- design_overview(config)
  bytes <- sum(2 * design$n * design$m * design$repetitions * 8)
  data.frame(datasets = nrow(design) * config$repetitions,
             method_fits = nrow(design) * config$repetitions * length(config$methods),
             raw_matrix_GiB = round(bytes / 1024^3, 3))
}
