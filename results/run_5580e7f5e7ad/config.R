list(schema_version = 1L, preset = "source_demo", seed = 20261004L, 
    m = 1000L, n = 300L, tau = 1.9, repetitions = 1L, scenarios = c("linear", 
    "confounded", "binary"), mixtures = list(sparse = c(H00 = 0.88, 
    H10 = 0.05, H01 = 0.05, H11 = 0.02), dense = c(H00 = 0.4, 
    H10 = 0.2, H01 = 0.2, H11 = 0.2)), dgp_profile = "source_code", 
    exposure_probability = 0.1, alpha_mean_scale = 0.05, beta_mean_scale = -0.5, 
    alpha_noise_variance = 1, beta_noise_variance = 4, direct_mean = 0.5, 
    direct_sd = 1, error_sd_m = 1, error_sd_y = 1, confounder_max = 0.5, 
    q = 0.05, methods = c("HDMT", "MDACT", "MLFDR"), method_options = list(
        HDMT = list(exact = 0L), MDACT = list(), MLFDR = list(
            eps = 0.01, twostep = FALSE, verbose = FALSE)), workers = 2L, 
    export_first_long_csv = FALSE)
