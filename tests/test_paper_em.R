PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
set.seed(20261015)
m <- 3000L
n <- 100L
states <- sample(1:4, m, TRUE, c(.4, .2, .2, .2))
alpha <- rnorm(m, .7, .06) * (states %in% c(2, 4))
beta <- rnorm(m, -.6, .10) * (states %in% c(3, 4))
va <- runif(m, .02, .025)
vb <- runif(m, .005, .006)
a <- alpha + rnorm(m, sd = sqrt(va))
b <- beta + rnorm(m, sd = sqrt(vb))
fit <- fit_mlfdr_paper(a, b, va, vb, n, list(eps = 1e-6, n_starts = 2L))
stopifnot(fit$converged, all(diff(fit$history) >= -1e-6),
          all(fit$lfdr >= 0 & fit$lfdr <= 1), max(abs(rowSums(fit$z) - 1)) < 1e-12,
          max(abs(fit$pi - c(.4, .2, .2, .2))) < .05,
          abs(fit$mu - .7) < .05, abs(fit$theta + .6) < .03)
raw <- fit_mlfdr_paper(a, b, va, vb, n,
                       list(eps = 1e-6, n_starts = 2L, scale = "raw", variance_tolerance = 1e-8))
stopifnot(max(abs(fit$lfdr - raw$lfdr)) < 1e-4,
          identical(lfdr_step_up(fit$lfdr, .05), lfdr_step_up(raw$lfdr, .05)))
e <- data.frame(alpha_hat = a, beta_hat = b, var_alpha = va, var_beta = vb)
result <- run_mlfdr(list(estimates = e, meta = list(n = n)), .05,
                    list(engine = "paper_em", eps = 1e-6, n_starts = 2L))
stopifnot(isTRUE(all.equal(result$score, fit$lfdr)),
          max(result$score[result$reject]) > .05,
          mean(result$score[result$reject]) <= .05)
cat("PASS: paper mixture recovery, monotone likelihood, valid posterior, scale invariance, adaptive step-up.\n")
