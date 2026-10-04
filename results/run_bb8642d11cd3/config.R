list(schema_version = 1L, preset = "teaching", seed = 20261010L, 
    m = 1000L, n = 300L, tau = c(1.3, 1.5, 1.9), repetitions = 40L, 
    scenarios = c("linear", "confounded", "binary"), mixtures = list(
        sparse = c(H00 = 0.88, H10 = 0.05, H01 = 0.05, H11 = 0.02
        ), dense = c(H00 = 0.4, H10 = 0.2, H01 = 0.2, H11 = 0.2
        )), dgp_profile = "teaching_fixed", exposure_probability = 0.5, 
    alpha_mean_scale = 0.2, beta_mean_scale = 0.5, alpha_noise_variance = 0, 
    beta_noise_variance = 0, direct_mean = 0.2, direct_sd = 0, 
    error_sd_m = 1, error_sd_y = 1, confounder_max = 0.3, q = 0.05, 
    methods = "ASHMED_Adaptive", method_options = list(ASHMED_Adaptive = list(
        grid_mode = "same", null_pseudocount = 0)), workers = 4L, 
    export_first_long_csv = TRUE, comparison_stage = "structure", 
    analysed_replicates = 1:3)
