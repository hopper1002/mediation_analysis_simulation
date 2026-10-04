PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
config <- default_config("teaching")
validate_config(config)
original <- default_config("paper_comparison")
for (key in c("n", "m", "tau", "repetitions", "scenarios", "mixtures", "q", "methods", "method_options"))
  stopifnot(identical(config[[key]], original[[key]]))
set.seed(123)
old_rng <- .Random.seed
grid <- build_design(config)
for (scenario in config$scenarios) for (mix in names(config$mixtures)) {
  idx <- which(grid$scenario == scenario & grid$mixture == mix)[1]
  data <- simulate_dataset(grid[idx, ], config)
  stopifnot(identical(old_rng, .Random.seed),
            identical(data, simulate_dataset(grid[idx, ], config)),
            sum(data$X == 1L) == 150L, sum(data$X == 0L) == 150L,
            all(table(factor(data$truth$state, levels=names(config$mixtures[[mix]]))) == config$m * config$mixtures[[mix]]),
            all(data$truth$alpha[data$truth$state %in% c("H00", "H01")] == 0),
            all(data$truth$beta[data$truth$state %in% c("H00", "H10")] == 0),
            all(data$truth$alpha[data$truth$state %in% c("H10", "H11")] == .2 * data$meta$tau),
            all(data$truth$beta[data$truth$state %in% c("H01", "H11")] == .5 * data$meta$tau))
  estimates <- estimate_coefficients(data)$estimates
  input <- make_method_input(data, estimates)
  stopifnot(!("truth" %in% names(input)), !anyNA(estimates))
  j <- which(data$truth$state == "H11")[1]
  df <- data.frame(M=data$M[,j], Y=data$Y[,j], X=data$X)
  if (length(data$Z)) df$Z <- data$Z
  fM <- if (length(data$Z)) M ~ X + Z else M ~ X
  fY <- if (length(data$Z)) Y ~ M + X + Z else Y ~ M + X
  a <- coef(summary(lm(fM, df)))["X", ]
  b <- if (scenario == "binary") coef(summary(glm(fY, df, family=binomial())))["M", ] else
    coef(summary(lm(fY, df)))["M", ]
  stopifnot(abs(estimates$alpha_hat[j] - a[1]) < 1e-10,
            abs(estimates$var_alpha[j] - a[2]^2) < 1e-10,
            abs(estimates$beta_hat[j] - b[1]) < 1e-10,
            abs(estimates$var_beta[j] - b[2]^2) < 1e-10)
}
rounded <- allocate_fixed_counts(101L, config$mixtures$sparse)
stopifnot(sum(rounded) == 101L, all(abs(rounded - 101 * config$mixtures$sparse) < 1))
invalid <- config
invalid$direct_sd <- 0.1
stopifnot(inherits(try(validate_config(invalid), silent=TRUE), "try-error"))
invalid <- config
invalid$alpha_mean_scale <- 0
stopifnot(inherits(try(validate_config(invalid), silent=TRUE), "try-error"))
cat("PASS: teaching DGP, fixed counts, RNG isolation, unchanged comparison settings, no truth input, OLS/GLM agreement.\n")
