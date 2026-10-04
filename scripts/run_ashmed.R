# Run the third comparison from already saved data; no data generation.
PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source(file.path(PROJECT_ROOT, "R/bootstrap.R"), encoding = "UTF-8")
source(file.path(PROJECT_ROOT, "R/ashmed_comparison.R"), encoding = "UTF-8")
opts <- list(ASHMED = list(J = 12L, r_min = .25, r_max = 8, correlations = c(-.5, 0, .5),
            initial_pi = c(H00 = .85, H10 = .1, H01 = .1, H11 = .05),
            relative_tolerance = 1e-8, kkt_tolerance = 1e-4, max_iter = 5000L,
            accelerate = TRUE, fit_subset = NULL, posterior_chunk_size = 5000L))
inputs <- read_ashmed_inputs(PROJECT_ROOT, method_options = opts)
suite <- run_ashmed_suite(inputs, PROJECT_ROOT)
for (label in names(suite)) {
  cat(label, ": ", suite[[label]]$experiment$output_dir, "\n", sep = "")
  print(suite[[label]]$fits[!suite[[label]]$fits$converged, ])
}
