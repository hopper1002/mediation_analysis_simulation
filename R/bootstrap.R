# Usage: source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
if (!exists("PROJECT_ROOT", inherits = TRUE)) stop("Set PROJECT_ROOT before sourcing bootstrap.R")
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
for (file in c("utils.R", "config.R", "simulation.R", "estimation.R", "methods/registry.R", "methods/mlfdr_em.R", "methods/builtin.R", "experiment.R", "reporting.R")) {
  sys.source(file.path(PROJECT_ROOT, "R", file), envir = environment())
}
mdact_vendor <- new.env(parent = asNamespace("stats"))
sys.source(file.path(PROJECT_ROOT, "vendor", "mdact_core.R"), envir = mdact_vendor)
