# Verify persisted decisions against raw truth and regenerate a delivery report.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
PROJECT_ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), winslash = "/")
source(file.path(PROJECT_ROOT, "R", "bootstrap.R"), encoding = "UTF-8")
directories <- list.dirs(file.path(PROJECT_ROOT, "results"), recursive = FALSE, full.names = TRUE)
directories <- directories[file.exists(file.path(directories, "experiment.rds"))]
matches <- Filter(function(path) {
  config <- readRDS(file.path(path, "config.rds"))
  config$dgp_profile == "paper2026" && config$m == 1000L && config$repetitions == 5L &&
    identical(config$methods, c("HDMT", "MDACT", "MLFDR"))
}, directories)
if (!length(matches)) stop("No completed demo results found.")
main_dir <- matches[[which.max(vapply(matches, function(p) as.numeric(file.info(file.path(p, "experiment.rds"))$mtime), numeric(1)))]]
saved <- readRDS(file.path(main_dir, "experiment.rds"))
config <- readRDS(file.path(main_dir, "config.rds"))
provenance <- readRDS(file.path(main_dir, "provenance.rds"))
manifest <- read_manifest(file.path(PROJECT_ROOT, "data", saved$data_id, "manifest.rds"))
metrics <- saved$metrics
stopifnot(nrow(manifest$design) == 180L, nrow(metrics) == 540L,
          all(metrics$status == "ok"), !anyNA(metrics$FDP), !anyNA(metrics$power),
          all(metrics$FDP >= 0 & metrics$FDP <= 1), all(metrics$power >= 0 & metrics$power <= 1))
stopifnot(identical(provenance$code, code_signature(PROJECT_ROOT)))
rounded <- 0L
for (id in manifest$design$case_id) {
  raw <- read_dataset(manifest, id)
  rows <- metrics[metrics$case_id == id, ]
  for (i in seq_len(nrow(rows))) {
    outcome <- readRDS(file.path(PROJECT_ROOT, "results", "method_cache", rows$method[i], paste0(rows$cache_key[i], ".rds")))
    validate_method_result(outcome$result, config$m)
    m <- compute_metrics(outcome$result$reject, raw$truth)
    stopifnot(isTRUE(all.equal(as.numeric(rows$FDP[i]), m$FDP)),
              isTRUE(all.equal(as.numeric(rows$power[i]), m$power)),
              rows$discoveries[i] == m$discoveries, rows$false_discoveries[i] == m$false_discoveries)
    if (rows$method[i] == "MLFDR") rounded <- rounded + outcome$result$diagnostics$roundoff_clipped
  }
}
stopifnot(isTRUE(all.equal(saved$summary, summarise_metrics(metrics), check.attributes = TRUE)))
figures <- list.files(file.path(main_dir, "figures"), "[.]png$", full.names = TRUE)
stopifnot(length(figures) == 6L, all(file.info(figures)$size > 10000L))
cat("PASS: all 180 file checksums and 540 decisions/metrics, summaries, code provenance, six figures.\n")

source_matches <- Filter(function(path) {
  c <- readRDS(file.path(path, "config.rds"))
  c$dgp_profile == "source_code" && c$m == 1000L && c$repetitions == 1L
}, directories)
stopifnot(length(source_matches) > 0L)
source_dir <- source_matches[[which.max(vapply(source_matches, function(p) as.numeric(file.info(file.path(p, "experiment.rds"))$mtime), numeric(1)))]]
source <- readRDS(file.path(source_dir, "experiment.rds"))
stopifnot(nrow(source$metrics) == 18L, all(source$metrics$status == "ok"))

relative <- function(path) substring(path, nchar(PROJECT_ROOT) + 2L)
lines <- c(
  "# \u8fd0\u884c\u4e0e\u4ea4\u4ed8\u68c0\u67e5\u8bb0\u5f55", "",
  "\u6267\u884c\u65e5\u671f\uff1a2026-10-04\uff08Asia/Shanghai\uff09\u3002", "",
  paste0("\u4e3b\u5b9e\u9a8c\uff1a`", relative(main_dir), "`\u3002"),
  paste0("\u6570\u636e\u5165\u53e3\uff1a`data/", saved$data_id, "/manifest.rds`\u3002"),
  paste0("\u6e90\u7801\u8bbe\u5b9a\u5bf9\u7167\uff1a`", relative(source_dir), "`\u3002"), "",
  "## \u5b9e\u9645\u6267\u884c", "",
  "- \u4e3b\u5b9e\u9a8c180\u4efd\u6570\u636e\u3001540\u6b21\u65b9\u6cd5\u8c03\u7528\uff0c\u5168\u90e8\u6210\u529f\uff1b3\u573a\u666f\u00d72\u6df7\u5408\u6bd4\u4f8b\u00d72\u6837\u672c\u91cf\u00d73\u4e2a\u03c4\u00d75\u91cd\u590d\u3002",
  "- \u989d\u5916\u6e90\u7801\u8bbe\u5b9a\u5bf9\u71676\u4efd\u6570\u636e\u300118\u6b21\u65b9\u6cd5\u8c03\u7528\uff0c\u5168\u90e8\u6210\u529f\uff0c\u6bcf\u6761\u4ef61\u6b21\u3002",
  paste0("- \u4e3b\u5b9e\u9a8c\u6709\u8b66\u544a\u7684\u8c03\u7528\uff1a", sum(nzchar(metrics$warning)), "\uff1b\u5747\u4fdd\u5b58\u5728 diagnostics.csv\u3002"),
  paste0("- MLFDR \u6d6e\u70b9\u622a\u754c\u6b21\u6570\uff08\u9010\u8def\u5f84\u7d2f\u8ba1\uff09\uff1a", rounded, "\uff1b\u53ea\u63a5\u53d71e-12\u5185\u7684\u820d\u5165\u8bef\u5dee\u3002"),
  "- \u5df2\u6838\u67e5\u5168\u90e8\u4e3b\u6570\u636eMD5\u3001\u6bcf\u65b9\u6cd5\u62d2\u7edd\u5411\u91cf\u4e0e\u771f\u503c\u8ba1\u7b97\u7684FDP/power\u3001CSV\u6c47\u603b\u4e0e6\u5f20PNG\u3002",
  "- notebook\u9010\u5355\u5143\u6267\u884c\uff1bHTML\u5bfc\u51fa\u548c\u56fe\u5f62\u5e03\u5c40\u53e6\u7531\u4ea4\u4ed8\u68c0\u67e5\u6838\u67e5\u3002", "",
  "## \u4e3b\u5b9e\u9a8c\u7ed3\u679c\u7684\u89e3\u91ca", "",
  sprintf("\u5404\u6761\u4ef6\u7ecf\u9a8cFDR\u5747\u503c\u8303\u56f4\uff1a%.4f\u2013%.4f\uff1bpower\u5747\u503c\u8303\u56f4\uff1a%.4f\u2013%.4f\u3002", min(saved$summary$FDR_mean), max(saved$summary$FDR_mean),
          min(saved$summary$power_mean), max(saved$summary$power_mean)),
  "\u53d1\u8868\u7248\u6b63\u6587\u7684\u03b1\u6548\u5e94\u5f88\u5f31\uff0c\u672c\u6b21\u6f14\u793a\u7684power\u4f4e\u3002\u6e90\u7801\u672a\u7f29\u653e\u7684\u7cfb\u6570\u566a\u58f0\u4e0e\u6b63\u6587\u4e0d\u540c\uff0c\u989d\u5916\u5bf9\u7167\u7684\u68c0\u51fa\u660e\u663e\u4e0d\u540c\u3002",
  "5\u6b21\u91cd\u590d\u7684\u6ce2\u52a8\u8f83\u5927\uff1b\u67d0\u4e2a\u6761\u4ef6\u7684\u5747\u503c\u8d85\u51fa\u76ee\u6807q\u4e0d\u7b49\u4e8e\u7a33\u5b9a\u8fdd\u80cc\u6216\u9a8c\u8bc1\u6e10\u8fd1\u4fdd\u8bc1\u3002\u5355\u6b21\u6e90\u7801\u5bf9\u7167\u66f4\u4e0d\u80fd\u4f5cFDR\u7ed3\u8bba\u3002", "",
  "| \u6761\u4ef6 | \u65b9\u6cd5 | \u5355\u6b21FDP | \u5355\u6b21power | \u53d1\u73b0\u6570 |",
  "|:--|:--|--:|--:|--:|"
)
for (i in seq_len(nrow(source$metrics))) {
  d <- source$metrics[i, ]
  lines <- c(lines, sprintf("| %s/%s, n=%d, \u03c4=%.1f | %s | %.4f | %.4f | %d |",
                            d$scenario, d$mixture, d$n, d$tau, d$method, d$FDP, d$power, d$discoveries))
}
lines <- c(lines, "", "## \u6b63\u786e\u6027\u4e0e\u5de5\u7a0b\u68c0\u67e5", "",
           "- run_tests.R\uff1a\u786e\u5b9a\u6027\u79cd\u5b50\u3001\u6269\u5c55\u7f51\u683c\u4e0d\u6270\u52a8\u539f\u573a\u666f\u79cd\u5b50\u3001OLS\u4e0elm\u7b49\u4ef7\u3001\u72b6\u6001/\u65b9\u5dee\u751f\u6210\u3001step-up\u96f6/\u5168/\u5e76\u5217/\u6d6e\u70b9\u8fb9\u754c\u3001\u6240\u6709\u9002\u914d\u5668\u3002",
           "- check_pipeline.R\uff1a\u6587\u4ef6\u5f80\u8fd4\u3001\u4e32\u884c/Windows PSOCK\u3001\u7f13\u5b58\u547d\u4e2d\u3001\u589e\u52a0\u63d2\u4ef6\u590d\u7528\u3001q\u53d8\u5316\u7f13\u5b58\u9694\u79bb\u3001checksum\u62d2\u7edd\u3001\u65b9\u6cd5\u5931\u8d25\u4e3aNA\u3002",
           "- \u4f7f\u7528\u5df2\u5b89\u88c5HDMT\u548cMLFDR\uff1bMDACT\u4f7f\u7528\u6709\u6765\u6e90\u548c\u8bb8\u53ef\u8bc1\u8bb0\u5f55\u7684\u914d\u5957\u4ee3\u7801\u3002",
           "- \u539fMLFDR-main\u672a\u4fee\u6539\uff1b\u9879\u76ee\u8fd0\u884c\u65e0\u9700\u8be5\u539f\u4ed3\u5e93\u8def\u5f84\u3002", "",
           "\u73af\u5883\u7248\u672c\u89c1\u672c\u6b21\u7ed3\u679c\u76ee\u5f55\u7684 package_versions.csv \u548c sessionInfo.txt\u3002\u6b64\u8bb0\u5f55\u662f\u672c\u673a\u8fd0\u884c\u9a8c\u8bc1\uff0c\u4e0d\u627f\u8bfa\u5176\u4ed6\u5305\u7248\u672c\u4f4d\u7ea7\u590d\u73b0\u3002")
writeLines(enc2utf8(lines), file.path(PROJECT_ROOT, "docs", "RUN_REPORT.md"), useBytes = TRUE)
writeLines(c(paste0("main_results=", relative(main_dir)), paste0("data=", relative(manifest$directory)),
             paste0("source_results=", relative(source_dir))), file.path(PROJECT_ROOT, "docs", "latest_run.txt"))
cat("Run report:", file.path(PROJECT_ROOT, "docs", "RUN_REPORT.md"), "\n")
