PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
set.seed(20261007)
for (i in seq_len(100)) {
  weights <- pmax(.001, stats::runif(3))
  weights <- weights / sum(weights)
  a <- weights[1]; b <- weights[2]; c <- weights[3]
  thresholds <- c(0, 10^runif(1, -9, -2), runif(1, 0, 1), 1)
  values <- vapply(sort(thresholds), mdact_vendor$mdact_null_cdf, numeric(1), a=a, b=b, c=c)
  stopifnot(all(values >= 0 & values <= 1), all(diff(values) >= -1e-12))
  for (t in thresholds[2:3]) {
    integrand <- function(v) {
      c0 <- c * v^2 + a * v
      c1 <- t - b * v
      s1 <- pmin(1, pmax(0, (t - c * v^2 - b * v) / a))
      s2 <- pmin(1, pmax(0, 2 * pmax(0, t - b * v) / (a + sqrt(a^2 + 4 * c * pmax(0, t - b * v)))))
      ifelse(c0 >= c1, s1, s2)
    }
    cut <- 2*t/(a+b+sqrt((a+b)^2+4*c*t))
    end <- 2*t/(b+sqrt(b^2+4*c*t))
    knots <- sort(unique(pmin(1, pmax(0, c(0, 1, cut, end, (t-c-a)/b)))))
    numerical <- sum(vapply(seq_len(length(knots)-1), function(j)
      integrate(integrand, knots[j], knots[j+1], subdivisions=2000,
                 abs.tol=1e-12, rel.tol=1e-10)$value, numeric(1)))
    analytical <- mdact_vendor$mdact_null_cdf(t, a, b, c)
    stopifnot(abs(numerical - analytical) < 1e-9)
  }
}
cat("PASS: MDACT analytic F00 agrees with independent numerical integration, including small tails.\n")
