# Base-R correctness checks: no testthat dependency.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash = "/")
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
check_dependencies()
config <- default_config()
config$n <- 120L
config$m <- 1000L
config$scenarios <- c("linear", "confounded", "binary")
config$tau <- 1.9
config$repetitions <- 1L
config$workers <- 1L
design <- build_design(config)
stopifnot(identical(build_design(config), design))
larger <- config
larger$repetitions <- 2L
larger$tau <- c(0.1, 1.9)
d2 <- build_design(larger)
for (i in seq_len(nrow(design))) {
  match_row <- d2$scenario == design$scenario[i] & d2$mixture == design$mixture[i] &
    d2$tau == design$tau[i] & d2$replicate == design$replicate[i]
  stopifnot(d2$seed[match_row] == design$seed[i])
}
cat("PASS: deterministic seeds remain stable when the grid expands.\n")

for (scenario in c("linear", "confounded")) {
  case <- design[which(design$scenario == scenario & design$mixture == "dense"), ]
  data <- simulate_dataset(case, config)
  stopifnot(identical(data, simulate_dataset(case, config)))
  fit <- estimate_coefficients(data)
  e <- fit$estimates
  for (j in c(1L, 11L, 100L)) {
    d <- data.frame(M = data$M[, j], Y = data$Y[, j], X = data$X)
    if (scenario == "confounded") d$Z <- data$Z
    fa <- if (scenario == "confounded") M ~ X + Z else M ~ X
    fb <- if (scenario == "confounded") Y ~ M + X + Z else Y ~ M + X
    a <- coef(summary(lm(fa, d)))["X", ]
    b <- coef(summary(lm(fb, d)))["M", ]
    stopifnot(isTRUE(all.equal(e$alpha_hat[j], unname(a[1]), tolerance = 1e-10)),
              isTRUE(all.equal(e$var_alpha[j], unname(a[2]^2), tolerance = 1e-10)),
              isTRUE(all.equal(e$p_alpha[j], max(unname(a[4]), 1e-17), tolerance = 1e-10)),
              isTRUE(all.equal(e$beta_hat[j], unname(b[1]), tolerance = 1e-10)),
              isTRUE(all.equal(e$var_beta[j], unname(b[2]^2), tolerance = 1e-10)),
              isTRUE(all.equal(e$p_beta[j], max(unname(b[4]), 1e-17), tolerance = 1e-10)))
  }
  input <- make_method_input(data, e)
  stopifnot(!("truth" %in% names(input)), !any(c("state", "alpha", "beta", "is_nonnull") %in% names(e)),
            !any(c("mixture", "tau", "seed", "profile") %in% names(input$meta)))
}
cat("PASS: vectorized OLS agrees with lm, including variances and p-values; no truth leakage.\n")

stopifnot(!any(lfdr_step_up(c(0.3, 0.4), 0.05)),
          all(lfdr_step_up(c(0.001, 0.002), 0.05)),
          identical(lfdr_step_up(c(0.01, 0.08, 0.9), 0.05), c(TRUE, TRUE, FALSE)),
          identical(lfdr_step_up(c(0.001, 0.001, 0.2, 0.2), 0.06), c(TRUE, TRUE, FALSE, FALSE)))
stopifnot(identical(lfdr_step_up(c(0.001, 1 + 2e-16), 0.05), c(TRUE, FALSE)),
          inherits(try(bound_lfdr(c(0.001, 1.1)), silent = TRUE), "try-error"))
stopifnot(inherits(try(validate_method_result(list(reject = c(1, 0)), 2), silent = TRUE), "try-error"))
truth <- data.frame(is_nonnull = c(TRUE, FALSE, TRUE, FALSE))
metrics <- compute_metrics(c(TRUE, TRUE, FALSE, FALSE), truth)
stopifnot(metrics$FDP == 0.5, metrics$power == 0.5,
          compute_metrics(rep(FALSE, 4), truth)$FDP == 0,
          is.na(compute_metrics(c(FALSE, FALSE), data.frame(is_nonnull = c(FALSE, FALSE)))$power))
cat("PASS: step-up edge cases, plugin output validation, FDP/power definitions.\n")

case <- design[which(design$scenario == "linear" & design$mixture == "dense"), ]
data <- simulate_dataset(case, config)
input <- make_method_input(data, estimate_coefficients(data)$estimates)
registry <- load_method_registry(PROJECT_ROOT)
stopifnot(all(c("HDMT", "MDACT", "MLFDR", "MaxP_BH") %in% names(registry)))
for (method in names(registry)) {
  captured <- capture_conditions(function() registry[[method]]$run(input, 0.05, config$method_options[[method]]))
  if (inherits(captured$value, "method_error")) stop(method, ": ", captured$value$message)
  validate_method_result(captured$value, config$m)
  cat(sprintf("PASS: %s adapter (%0.3f seconds).\n", method, captured$elapsed))
}

# Distribution check distinguishes variance/n from source sd=4 and confirms the true state masks.
big <- config
big$m <- 20000L
big$mixtures <- list(all_signal = c(H00 = 0, H10 = 0, H01 = 0, H11 = 1))
big$scenarios <- "linear"
big$n <- 100L
b <- simulate_dataset(build_design(big)[1, ], big)
stopifnot(abs(var(b$truth$alpha) - 1/100) < 0.0005,
          abs(var(b$truth$beta) - 4/100) < 0.002)
cat("PASS: paper-profile coefficient noise variances are 1/n and 4/n.\n")
cat("All correctness checks passed.\n")
