# Vectorized Frisch-Waugh-Lovell OLS, with an intercept in both regressions.
# Algebraically identical to lm(M_i ~ X [+ Z]) and lm(Y_i ~ M_i + X [+ Z]).
estimate_coefficients <- function(data) {
  n <- nrow(data$M)
  m <- ncol(data$M)
  C <- cbind(intercept = 1, X = data$X)
  if (length(data$Z)) C <- cbind(C, Z = data$Z)
  qrC <- qr(C)
  if (qrC$rank != ncol(C)) stop("Rank-deficient exposure/covariate design.")
  bM <- qr.coef(qrC, data$M)
  residual_M <- qr.resid(qrC, data$M)
  df_alpha <- n - ncol(C)
  sigma2_M <- colSums(residual_M^2) / df_alpha
  invCtC <- chol2inv(chol(crossprod(C)))
  dimnames(invCtC) <- list(colnames(C), colnames(C))
  alpha_hat <- as.numeric(bM["X", ])
  var_alpha <- sigma2_M * invCtC["X", "X"]
  p_alpha <- 2 * stats::pt(-abs(alpha_hat / sqrt(var_alpha)), df = df_alpha)
  warnings <- character()
  if (data$meta$scenario != "binary") {
    residual_Y <- qr.resid(qrC, data$Y)
    ssM <- colSums(residual_M^2)
    beta_hat <- colSums(residual_M * residual_Y) / ssM
    residual_full <- residual_Y - sweep(residual_M, 2, beta_hat, "*")
    df_beta <- n - ncol(C) - 1L
    var_beta <- colSums(residual_full^2) / df_beta / ssM
    p_beta <- 2 * stats::pt(-abs(beta_hat / sqrt(var_beta)), df = df_beta)
  } else {
    beta_hat <- var_beta <- p_beta <- rep(NA_real_, m)
    for (j in seq_len(m)) {
      captured <- capture_conditions(function() {
        d <- data.frame(y = data$Y[, j], mediator = data$M[, j], X = data$X)
        fit <- stats::glm(y ~ mediator + X, data = d, family = stats::binomial())
        if (!fit$converged) stop("Logistic regression did not converge for pathway ", j)
        coef <- stats::coef(summary(fit))["mediator", ]
        c(beta = coef[[1]], variance = coef[[2]]^2, p = coef[[4]])
      })
      warnings <- c(warnings, captured$warnings)
      if (inherits(captured$value, "method_error")) stop(captured$value$message)
      beta_hat[j] <- captured$value[[1]]
      var_beta[j] <- captured$value[[2]]
      p_beta[j] <- captured$value[[3]]
    }
  }
  estimates <- data.frame(pathway_id = colnames(data$M), alpha_hat = alpha_hat, beta_hat = beta_hat,
                          var_alpha = var_alpha, var_beta = var_beta,
                          se_alpha = sqrt(var_alpha), se_beta = sqrt(var_beta),
                          p_alpha = pmax(p_alpha, 1e-17), p_beta = pmax(p_beta, 1e-17))
  numeric_columns <- estimates[, -1, drop = FALSE]
  if (!all(is.finite(as.matrix(numeric_columns))) || any(var_alpha <= 0 | var_beta <= 0))
    stop("Invalid coefficient estimates/standard errors; dataset not analysed.")
  list(estimates = estimates, warnings = unique(warnings))
}

make_method_input <- function(data, estimates) {
  # Ground truth is deliberately excluded from the method interface.
  list(estimates = estimates,
       raw = list(X = data$X, Z = data$Z, M = data$M, Y = data$Y),
       meta = data$meta[c("case_id", "n", "m", "scenario")])
}
