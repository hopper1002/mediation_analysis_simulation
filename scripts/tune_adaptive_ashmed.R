PROJECT_ROOT <- normalizePath(getwd(), winslash = "/")
source("R/bootstrap.R", encoding = "UTF-8")
source("R/ashmed_optimization.R", encoding = "UTF-8")
inputs <- optimization_inputs(PROJECT_ROOT)
dir <- optimization_directory(PROJECT_ROOT)
# Prespecified penalty strengths on the identifiable, shifted-only dictionary.
arms <- lapply(c(0, 10, 30, 90), function(lambda) list(zero_anchor = FALSE, null_pseudocount = lambda))
names(arms) <- paste0("null", c(0, 10, 30, 90))
atomic_save_rds(arms, file.path(dir, "regularization_options.rds"))
dput(arms, file = file.path(dir, "regularization_options.R"))
scores <- list()
for (arm in names(arms)) {
  message("Regularization development: ", arm)
  result <- run_adaptive_trial(inputs, PROJECT_ROOT, arms[[arm]], "development", 4:8, arm)
  scores[[arm]] <- score_adaptive_trial(result, arm)
  print(scores[[arm]], row.names = FALSE)
  write_csv(do.call(rbind, scores), file.path(dir, "development_scores.csv"))
}
table <- do.call(rbind, scores)
passing <- names(arms)[vapply(names(arms), function(a) all(table$passes[table$arm == a]), logical(1))]
if (!length(passing)) stop("No candidate passed the prespecified development gates; continue development before using held-out repeats.")
quality <- vapply(passing, function(a) mean(table$mean_power[table$arm == a]), numeric(1))
selected <- passing[which.max(quality)]
selection <- list(arm = selected, options = arms[[selected]], development_replicates = 1:8,
  confirmation_replicates = 9:12, validation_replicates = 13:40, scores = table,
  fingerprint = method_fingerprint(load_method_registry(PROJECT_ROOT)$ASHMED_Adaptive, PROJECT_ROOT))
atomic_save_rds(selection, file.path(dir, "selected_candidate.rds"))
dput(selection[c("arm", "options", "development_replicates", "confirmation_replicates", "validation_replicates")],
     file = file.path(dir, "selected_candidate.R"))
cat("Selected candidate: ", selected, "\n", sep = "")
confirm <- run_adaptive_trial(inputs, PROJECT_ROOT, selection$options, "confirmation", 9:12, selected)
score <- score_adaptive_trial(confirm, selected)
write_csv(score, file.path(dir, "confirmation_scores.csv"))
print(score, row.names = FALSE)
if (!all(score$passes)) stop("Confirmation failed. Do not use the validation subset for tuning.")
atomic_save_rds(list(selection = selection, confirmation = score), file.path(dir, "locked_method.rds"))
cat("Method locked; repeats 13–40 remain reserved for final evaluation.\n")
