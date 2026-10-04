list(schema_version = 1L, preset = "demo", seed = 20261004L, 
    m = 1000L, n = c(100L, 300L), tau = c(0.1, 1, 1.9), repetitions = 5L, 
    scenarios = c("linear", "confounded", "binary"), mixtures = list(
        sparse = c(H00 = 0.88, H10 = 0.05, H01 = 0.05, H11 = 0.02
        ), dense = c(H00 = 0.4, H10 = 0.2, H01 = 0.2, H11 = 0.2
        )), dgp_profile = "paper2026", exposure_probability = 0.1, 
    alpha_mean_scale = 0.05, beta_mean_scale = -0.5, alpha_noise_variance = 1, 
    beta_noise_variance = 4, direct_mean = 1, direct_sd = 0.707106781186548, 
    error_sd_m = 1, error_sd_y = 1, confounder_max = 0.5, q = 0.05, 
    methods = c("HDMT", "MDACT", "MLFDR"), method_options = list(
        HDMT = list(exact = 0L), MDACT = list(), MLFDR = list(
            eps = 0.01, twostep = FALSE, verbose = FALSE)), workers = 2L, 
    export_first_long_csv = TRUE)
