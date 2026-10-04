# Copy to R/methods/plugins/my_method.R, implement run(), then add "MyMethod"
# to config$methods in the notebook. Never read simulation truth inside a method.
register_plugin <- function(registry) {
  run <- function(input, q, options = list()) {
    estimates <- input$estimates
    # input$raw contains X, Z, M, Y for methods requiring original observations.
    # Replace this stop() with your algorithm, e.g. yourPackage::fit(...).
    stop("Implement MyMethod before activating it.")
    # Return one decision per pathway in the original input order:
    # list(reject = logical(nrow(estimates)), score = numeric(nrow(estimates)),
    #      score_type = "p_value", diagnostics = list(...))
  }
  register_method(registry, "MyMethod", run, packages = character(), version = "0.1.0",
                  description = "Replace with your method and citation")
}
