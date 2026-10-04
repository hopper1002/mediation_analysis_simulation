summarise_metrics <- function(metrics) {
  keys <- c("scenario", "mixture", "profile", "n", "m", "tau", "method", "q")
  groups <- split(seq_len(nrow(metrics)), interaction(metrics[, keys], drop = TRUE))
  rows <- lapply(groups, function(idx) {
    d <- metrics[idx, ]
    row <- d[1, keys, drop = FALSE]
    row$n_total <- nrow(d)
    row$n_success <- sum(d$status == "ok")
    row$n_error <- sum(d$status == "error")
    row$n_warning <- sum(nzchar(d$warning))
    # Empirical FDR = mean of per-replicate FDP, not pooled V / pooled R.
    for (metric in c("FDP", "power", "discoveries", "elapsed_seconds")) {
      x <- d[[metric]][d$status == "ok" & is.finite(d[[metric]])]
      prefix <- if (metric == "FDP") "FDR" else metric
      count <- length(x)
      mu <- if (count) mean(x) else NA_real_
      sd_x <- if (count > 1L) stats::sd(x) else NA_real_
      se <- sd_x / sqrt(count)
      half <- if (count > 1L) stats::qt(0.975, count - 1L) * se else NA_real_
      row[[paste0(prefix, "_n")]] <- count
      row[[paste0(prefix, "_mean")]] <- mu
      row[[paste0(prefix, "_sd")]] <- sd_x
      row[[paste0(prefix, "_mcse")]] <- se
      row[[paste0(prefix, "_lo")]] <- if (prefix %in% c("FDR", "power")) max(0, mu - half) else mu - half
      row[[paste0(prefix, "_hi")]] <- if (prefix %in% c("FDR", "power")) min(1, mu + half) else mu + half
      # All-zero/all-one bounded samples do not justify a zero-width interval.
      # A worst-case one-sided endpoint bound, with 0.025 per tail, is valid
      # for any iid observations in [0,1], including per-replicate FDP/power.
      if (prefix %in% c("FDR", "power") && count > 1L && all(x == 0))
        row[[paste0(prefix, "_hi")]] <- 1 - 0.025^(1 / count)
      if (prefix %in% c("FDR", "power") && count > 1L && all(x == 1))
        row[[paste0(prefix, "_lo")]] <- 0.025^(1 / count)
    }
    row
  })
  result <- do.call(rbind, rows)
  # Remove column names so the "method" column is not bound to order(method=...).
  result <- result[do.call(order, unname(as.list(result[, keys]))), ]
  rownames(result) <- NULL
  result
}

plot_comparison <- function(summary, scenario, metric = c("FDR", "power")) {
  metric <- match.arg(metric)
  d <- summary[summary$scenario == scenario, , drop = FALSE]
  d$estimate <- d[[paste0(metric, "_mean")]]
  d$lower <- d[[paste0(metric, "_lo")]]
  d$upper <- d[[paste0(metric, "_hi")]]
  d$mixture <- factor(d$mixture, levels = c("sparse", "dense", setdiff(unique(d$mixture), c("sparse", "dense"))))
  colors <- c(HDMT = "#2864A0", MDACT = "#13856B", MLFDR = "#C1486A", MaxP_BH = "#9C6900")
  extra <- setdiff(unique(d$method), names(colors))
  if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  y_limits <- if (metric == "FDR") c(0, max(0.12, d$upper, na.rm = TRUE)) else
    c(0, min(1, 1.05 * max(0.05, d$upper, na.rm = TRUE)))
  graph <- ggplot2::ggplot(d, ggplot2::aes(x = tau, y = estimate, color = method, group = method)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = method), alpha = 0.1, color = NA) +
    ggplot2::geom_line(linewidth = 0.75) + ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(mixture ~ n, labeller = ggplot2::label_both) +
    ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$tau))) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(title = paste(scenario, "|", if (metric == "FDR") "Empirical FDR" else "Power"),
                  subtitle = paste0("m = ", unique(d$m), "; ", unique(d$profile),
                                    "; ribbons: approximate 95% Monte Carlo intervals"),
                  x = "Signal strength (tau)", y = if (metric == "FDR") "Mean FDP" else "Mean power",
                  color = "Method", fill = "Method") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"))
  if (metric == "FDR") graph <- graph + ggplot2::geom_hline(yintercept = unique(d$q), color = "#444444", linetype = 2)
  graph
}

save_figures <- function(experiment) {
  directory <- ensure_dir(file.path(experiment$output_dir, "figures"))
  paths <- character()
  for (scenario in unique(experiment$summary$scenario)) for (metric in c("FDR", "power")) {
    graph <- plot_comparison(experiment$summary, scenario, metric)
    path <- file.path(directory, paste0(scenario, "_", tolower(metric), ".png"))
    ggplot2::ggsave(path, graph, width = 9, height = 6, dpi = 160, bg = "white")
    paths <- c(paths, path)
  }
  invisible(paths)
}

summarise_paired_metrics <- function(metrics, left = "MLFDR", right = c("HDMT", "MDACT")) {
  keys <- c("scenario", "mixture", "profile", "n", "m", "tau", "q")
  output <- list()
  for (other in right) {
    lhs <- metrics[metrics$method == left & metrics$status == "ok", c("case_id", keys, "FDP", "power")]
    rhs <- metrics[metrics$method == other & metrics$status == "ok", c("case_id", "FDP", "power")]
    pairs <- merge(lhs, rhs, by = "case_id", suffixes = c("_left", "_right"))
    groups <- split(seq_len(nrow(pairs)), interaction(pairs[, keys], drop = TRUE))
    for (idx in groups) {
      d <- pairs[idx, ]
      row <- d[1, keys, drop = FALSE]
      row$comparison <- paste(left, "-", other)
      row$n_pairs <- nrow(d)
      for (name in c("FDP", "power")) {
        delta <- d[[paste0(name, "_left")]] - d[[paste0(name, "_right")]]
        delta <- delta[is.finite(delta)]
        prefix <- if (name == "FDP") "delta_FDR" else "delta_power"
        se <- if (length(delta) > 1) sd(delta) / sqrt(length(delta)) else NA_real_
        half <- if (length(delta) > 1) qt(.975, length(delta) - 1) * se else NA_real_
        row[[paste0(prefix, "_mean")]] <- mean(delta)
        row[[paste0(prefix, "_mcse")]] <- se
        row[[paste0(prefix, "_lo")]] <- mean(delta) - half
        row[[paste0(prefix, "_hi")]] <- mean(delta) + half
      }
      output[[length(output) + 1L]] <- row
    }
  }
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}

plot_paper_comparison <- function(summary, scenario) {
  original <- summary[summary$scenario == scenario, ]
  d <- do.call(rbind, lapply(c("FDR", "power"), function(metric) {
    x <- original
    x$metric <- if (metric == "FDR") "FDR" else "Power"
    x$estimate <- x[[paste0(metric, "_mean")]]
    x$lower <- x[[paste0(metric, "_lo")]]
    x$upper <- x[[paste0(metric, "_hi")]]
    x
  }))
  conditions <- unique(original[, c("mixture", "n")])
  conditions <- conditions[order(match(conditions$mixture, c("sparse", "dense")), conditions$n), ]
  labels <- unlist(lapply(seq_len(nrow(conditions)), function(i)
    paste(conditions$mixture[i], paste0("n=", conditions$n[i]), c("FDR", "Power"), sep = " | ")))
  d$panel <- factor(paste(d$mixture, paste0("n=", d$n), d$metric, sep = " | "), levels = labels)
  anchors <- unique(d[, c("panel", "metric", "tau")])
  anchors$lo <- 0
  anchors$hi <- ifelse(anchors$metric == "Power", 1, max(.10, d$upper[d$metric == "FDR"], na.rm = TRUE))
  qline <- unique(d[d$metric == "FDR", c("panel", "q")])
  colors <- c(HDMT = "#2864A0", MDACT = "#13856B", MLFDR = "#C1486A", MaxP_BH = "#9C6900")
  extra <- setdiff(unique(d$method), names(colors))
  if (length(extra)) colors <- c(colors, setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  ggplot2::ggplot(d, ggplot2::aes(tau, estimate, color = method, group = method)) +
    ggplot2::geom_blank(data = anchors, ggplot2::aes(x = tau, y = lo), inherit.aes = FALSE) +
    ggplot2::geom_blank(data = anchors, ggplot2::aes(x = tau, y = hi), inherit.aes = FALSE) +
    ggplot2::geom_hline(data = qline, ggplot2::aes(yintercept = q), inherit.aes = FALSE,
                        color = "#555555", linetype = 2) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = method), alpha = .09, color = NA) +
    ggplot2::geom_line(linewidth = .85) + ggplot2::geom_point(size = 2.4) +
    ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(breaks = sort(unique(d$tau))) +
    ggplot2::labs(title = paste(scenario, "| FDR and power"),
                  subtitle = paste0("m=", unique(d$m), "; B=", paste(unique(d$n_success), collapse = ","),
                                    "; same saved datasets for all methods"),
                  x = "Signal strength (tau)", y = "Estimate", color = "Method", fill = "Method") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"))
}

collect_fit_diagnostics <- function(experiment, root) {
  rows <- experiment$metrics[experiment$metrics$method == "MLFDR" & experiment$metrics$status == "ok", ]
  output <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, ]
    outcome <- readRDS(file.path(root, "results", "method_cache", "MLFDR", paste0(row$cache_key, ".rds")))
    fit <- outcome$result$diagnostics
    if (is.null(fit$engine) || fit$engine != "paper_em") return(NULL)
    data.frame(case_id = row$case_id, converged = fit$converged, iterations = fit$iterations,
               mu = fit$mu, theta = fit$theta, kappa = fit$kappa, psi = fit$psi,
               pi00 = fit$pi[1], pi10 = fit$pi[2], pi01 = fit$pi[3], pi11 = fit$pi[4],
               selected_start = fit$selected_start, loglik = fit$loglik)
  })
  result <- do.call(rbind, output)
  if (!is.null(result)) rownames(result) <- NULL
  result
}

save_comparison_artifacts <- function(experiment, root) {
  directory <- ensure_dir(file.path(experiment$output_dir, "figures"))
  for (scenario in unique(experiment$summary$scenario))
    ggplot2::ggsave(file.path(directory, paste0(scenario, "_comparison.png")),
                    plot_paper_comparison(experiment$summary, scenario),
                    width = 10, height = 7, dpi = 160, bg = "white")
  paired <- summarise_paired_metrics(experiment$metrics)
  write_csv(paired, file.path(experiment$output_dir, "paired_differences.csv"))
  fits <- collect_fit_diagnostics(experiment, root)
  if (!is.null(fits)) write_csv(fits, file.path(experiment$output_dir, "fit_diagnostics.csv"))
  invisible(list(paired = paired, fits = fits))
}
