register_method <- function(registry, name, fn, packages = character(), version = "1.0.0", description = "") {
  stopifnot(is.list(registry), is.function(fn), length(name) == 1L,
            grepl("^[A-Za-z][A-Za-z0-9_]*$", name), !(name %in% names(registry)))
  registry[[name]] <- list(name = name, run = fn, packages = packages,
                           version = version, description = description)
  registry
}

validate_method_result <- function(result, m) {
  if (!is.list(result) || !is.logical(result$reject) || length(result$reject) != m || anyNA(result$reject))
    stop("Method must return list(reject = logical(m)), without NA.")
  if (!is.null(result$score) && (length(result$score) != m || !is.numeric(result$score) || any(!is.finite(result$score))))
    stop("Optional score must be a finite numeric vector of length m.")
  result
}

load_method_registry <- function(root) {
  registry <- list()
  registry <- register_method(registry, "HDMT", run_hdmt, "HDMT", "1.0.0", "HDMT::fdr_est (paper uses exact=0)")
  registry <- register_method(registry, "MDACT", run_mdact, "HDMT", "1.1.0", "Author MDACT FDR threshold with equivalent analytic null CDF")
  registry <- register_method(registry, "MLFDR", run_mlfdr, "MLFDR", "2.0.0", "Paper four-state EM or MLFDR::localFDR + adaptive step-up; engine in config")
  plugins <- sort(list.files(file.path(root, "R", "methods", "plugins"), "[.]R$", full.names = TRUE))
  for (path in plugins) {
    env <- new.env(parent = environment())
    sys.source(path, envir = env)
    if (!exists("register_plugin", envir = env, inherits = FALSE)) stop("Plugin missing register_plugin(): ", path)
    old_names <- names(registry)
    registry <- env$register_plugin(registry)
    for (name in setdiff(names(registry), old_names)) registry[[name]]$source_hash <- file_hash(path)
  }
  registry
}

registry_info <- function(registry) {
  do.call(rbind, lapply(registry, function(x) data.frame(
    method = x$name, adapter_version = x$version, packages = paste(x$packages, collapse = ", "),
    description = x$description, stringsAsFactors = FALSE)))
}
