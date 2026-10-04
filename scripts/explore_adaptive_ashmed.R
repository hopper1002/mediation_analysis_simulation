PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source("R/bootstrap.R", encoding = "UTF-8")
source("R/ashmed_optimization.R", encoding = "UTF-8")
inputs <- optimization_inputs(PROJECT_ROOT)
arms <- list(
  zero_rect = list(learn_location = FALSE, zero_anchor = TRUE, J = 8L, r_max = 16, correlations = 0, null_pseudocount = 0),
  shifted_same = list(grid_mode = "same", zero_anchor = TRUE, null_pseudocount = 0),
  shifted_rect = list(zero_anchor = TRUE, null_pseudocount = 0),
  shifted_rect_no_anchor = list(zero_anchor = FALSE, null_pseudocount = 0)
)
dir <- optimization_directory(PROJECT_ROOT)
atomic_save_rds(arms, file.path(dir, "structure_options.rds"))
dput(arms, file = file.path(dir, "structure_options.R"))
scores <- list()
for (arm in names(arms)) {
  message("Structure exploration: ", arm)
  result <- run_adaptive_trial(inputs, PROJECT_ROOT, arms[[arm]], "structure", 1:3, arm)
  scores[[arm]] <- score_adaptive_trial(result, arm)
  print(scores[[arm]], row.names = FALSE)
  write_csv(do.call(rbind, scores), file.path(dir, "structure_scores.csv"))
}
