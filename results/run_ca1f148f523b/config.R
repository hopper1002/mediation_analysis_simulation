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
        J = 6L, r_min = 0.25, r_max = 8, correlations = 0, grid_mode = "rectangular", 
        learn_location = TRUE, zero_anchor = FALSE, include_atom = TRUE, 
        minimum_location_se = 0.5, null_pseudocount = 10, relative_tolerance = 1e-08, 
        kkt_tolerance = 1e-04, max_iter = 5000L, accelerate = TRUE, 
        location_max_iter = 1500L, location_tolerance = 1e-06)), 
    workers = 4L, export_first_long_csv = TRUE, comparison_stage = "development", 
    analysed_replicates = 4:8)
