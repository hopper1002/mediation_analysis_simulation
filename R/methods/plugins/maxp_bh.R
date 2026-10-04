# Runnable extension example; inactive unless config$methods includes "MaxP_BH".
register_plugin <- function(registry) {
  run <- function(input, q, options = list()) {
    e <- input$estimates
    adjusted <- stats::p.adjust(pmax(e$p_alpha, e$p_beta), method = "BH")
    list(reject = adjusted <= q, score = adjusted, score_type = "BH_adjusted_p",
         diagnostics = list(rule = "Max-P followed by Benjamini-Hochberg"))
  }
  register_method(registry, "MaxP_BH", run, packages = "stats", version = "1.0.0",
                  description = "Runnable plugin example: Max-P + BH")
}
