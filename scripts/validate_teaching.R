PROJECT_ROOT <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
source(file.path(PROJECT_ROOT,"R","bootstrap.R"),encoding="UTF-8")
pointer <- readRDS(file.path(PROJECT_ROOT,"docs","current_teaching.rds"))
config <- pointer$config
validate_config(config)
stopifnot(config$dgp_profile == "teaching_fixed")
manifest <- read_manifest(pointer$manifest_path)
design <- build_design(config)
stopifnot(identical(design$case_id,manifest$design$case_id),
          identical(design$seed,manifest$design$seed), identical(manifest$dgp,dgp_config(config)))
generator <- code_signature(PROJECT_ROOT,file.path(PROJECT_ROOT,"R",c("utils.R","config.R","simulation.R")))
expected <- paste0("data_",substr(object_hash(list(dgp=dgp_config(config),generator=generator)),1,12))
stopifnot(identical(expected,manifest$data_id))
run_dir <- pointer$output_dir
record <- readRDS(file.path(run_dir,"experiment.rds"))
metrics <- record$metrics
stopifnot(nrow(metrics)==nrow(design)*length(config$methods),all(metrics$status=="ok"),
          all(metrics$q==config$q),!anyDuplicated(paste(metrics$case_id,metrics$method)),
          all(metrics$FDP>=0 & metrics$FDP<=1),all(metrics$power>=0 & metrics$power<=1))
for (i in seq_len(nrow(design))) {
  data <- read_dataset(manifest,design$case_id[i])
  pi <- config$mixtures[[design$mixture[i]]]
  actual_counts <- table(factor(data$truth$state,levels=names(pi)))
  stopifnot(all(actual_counts==allocate_fixed_counts(config$m,pi)),
            sum(data$X==1)==round(data$meta$n*config$exposure_probability),
            all(data$truth$alpha==config$alpha_mean_scale*data$meta$tau*(data$truth$state %in% c("H10","H11"))),
            all(data$truth$beta==config$beta_mean_scale*data$meta$tau*(data$truth$state %in% c("H01","H11"))),
            all(data$truth$direct==config$direct_mean))
  rows <- metrics[metrics$case_id==data$meta$case_id,]
  for (j in seq_len(nrow(rows))) {
    row <- rows[j,]
    outcome <- readRDS(file.path(PROJECT_ROOT,"results","method_cache",row$method,paste0(row$cache_key,".rds")))
    actual <- compute_metrics(outcome$result$reject,data$truth)
    for (name in names(actual))
      stopifnot(isTRUE(all.equal(as.numeric(row[[name]]),as.numeric(actual[[name]]),tolerance=1e-12)))
  }
}
stopifnot(isTRUE(all.equal(summarise_metrics(metrics),record$summary,tolerance=1e-12)))
paired <- utils::read.csv(file.path(run_dir,"paired_differences.csv"))
stopifnot(isTRUE(all.equal(summarise_paired_metrics(metrics),paired,tolerance=1e-12,check.attributes=FALSE)))
fits <- utils::read.csv(file.path(run_dir,"fit_diagnostics.csv"))
stopifnot(nrow(fits)==nrow(design),all(fits$converged))
provenance <- readRDS(file.path(run_dir,"provenance.rds"))
stopifnot(identical(provenance$code,code_signature(PROJECT_ROOT)))
for (scenario in config$scenarios)
  stopifnot(file.exists(file.path(run_dir,"figures",paste0(scenario,"_comparison.png"))))
probe <- analyse_case(1L,manifest,config,PROJECT_ROOT,load_method_registry(PROJECT_ROOT))
stopifnot(all(probe$status=="ok"),all(probe$cached))
writeLines(c(paste0("run_dir=",run_dir),paste0("manifest_path=",pointer$manifest_path)),
           file.path(PROJECT_ROOT,"docs","current_teaching.txt"))
keys <- c("seed","m","n","tau","repetitions","exposure_probability","alpha_mean_scale", "beta_mean_scale",
          "direct_mean","error_sd_m","error_sd_y","confounder_max","q")
metadata <- data.frame(parameter=keys,value=vapply(config[keys],function(x) paste(x,collapse="|"),character(1)))
write_csv(metadata,file.path(run_dir,"teaching_design.csv"))
cat(sprintf("PASS: %d saved datasets, fixed state counts/exposure, %d decisions and metrics, summaries, paired differences, EM convergence, provenance, figures, CLI/IRkernel parity.\n",nrow(design),nrow(metrics)))
cat("Run:",run_dir,"\n")
