# Project utilities. No global workspace clearing or working-directory changes.
ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

object_hash <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  # v3 embeds native_encoding in its header (C vs UTF-8 differs between
  # Windows Rscript and IRkernel). v2 makes ASCII/numeric seed and cache
  # configuration hashes independent of that locale metadata.
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}

atomic_save_rds <- function(x, path, compress = TRUE) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = "write-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = compress, version = 3)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Cannot save: ", path)
  invisible(path)
}

write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  invisible(path)
}

file_hash <- function(path) unname(tools::md5sum(path))

code_signature <- function(root, files = NULL) {
  if (is.null(files)) {
    files <- sort(c(list.files(file.path(root, "R"), "[.]R$", recursive = TRUE, full.names = TRUE),
                    list.files(file.path(root, "vendor"), "[.]R$", recursive = TRUE, full.names = TRUE)))
  }
  setNames(vapply(files, file_hash, character(1)), substring(files, nchar(root) + 2L))
}

with_seed <- function(seed, fn) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  fn()
}

capture_conditions <- function(fn) {
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(fn(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(list(message = conditionMessage(e)), class = "method_error"))
  list(value = value, warnings = unique(warnings), elapsed = proc.time()[["elapsed"]] - started)
}

package_versions <- function(packages) {
  data.frame(package = packages, version = vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1)), stringsAsFactors = FALSE)
}

check_dependencies <- function() {
  versions <- package_versions(c("HDMT", "MLFDR", "NMOF", "ggplot2", "IRkernel"))
  required <- versions$package != "IRkernel"
  if (any(is.na(versions$version[required])))
    stop("Missing R packages: ", paste(versions$package[required & is.na(versions$version)], collapse = ", "))
  versions
}
