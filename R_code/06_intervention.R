# =============================================================================
# 06_intervention.R
# Intervention dose-response: shock-only + shock+damping
# =============================================================================
message("=== 06_intervention.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
sector_names_intervention <- c("coal", "oil", "NGas", "thermal", "electric")
forecast_years  <- CFG$future_years
n_future        <- length(forecast_years)
intervention_year <- CFG$intervention_year
intervention_start_index <- which(forecast_years >= intervention_year)[1]
reduction_grid <- seq(0, 0.50, by = 0.05)
damping_grid   <- c(0, 0.25, 0.50, 0.75)
N_DRAWS <- 4000
set.seed(12345)
extract_forecast_parameters <- function(fit, n_draws = 4000) {
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y"))
  total <- length(post$s_mu)
  idx <- if (total > n_draws) sample(seq_len(total), n_draws, replace = FALSE) else seq_len(total)
  list(mu_T1 = post$mu_trend[idx, ncol(post$mu_trend) - 1], mu_T  = post$mu_trend[idx, ncol(post$mu_trend)],
       sigma_mu = post$s_mu[idx],sigma_Y  = post$s_Y[idx])
}

generate_common_randomness <- function(n_draws, n_future) {
  list(z_mu = matrix(rnorm(n_draws * n_future), n_draws, n_future),
       z_y  = matrix(rnorm(n_draws * n_future), n_draws, n_future))
}

simulate_intervention <- function(params, rnd, reduction_fraction, damping,start_idx, n_future) {
  n_draws <- length(params$sigma_mu)
  gamma <- log(1 - reduction_fraction)
  logY <- matrix(NA_real_, n_draws, n_future)
  prev <- params$mu_T1; curr <- params$mu_T
  for (h in seq_len(n_future)) {
    s_mu_h <- if (h < start_idx) params$sigma_mu else (1 - damping) * params$sigma_mu
    next_base <- 2 * curr - prev + s_mu_h * rnd$z_mu[, h]
    shift <- if (h >= start_idx) gamma else 0
    next_mu <- next_base + shift
    logY[, h] <- next_mu + params$sigma_Y * rnd$z_y[, h]
    prev <- curr; curr <- next_base
  }
  logY
}

step6_rows <- list()
counter <- 0
for (sector in sector_names_intervention) {
  pretty <- SECTOR_LABELS[[sector]]
  Y_obs  <- as.numeric(D[[sector]])
  Y_2023 <- tail(Y_obs, 1)
  target_log <- log(Y_2023) + log(CFG$target_ratio)
  
  fit <- step2_results$fits[[sector]]$Reduced_SORW
  params <- extract_forecast_parameters(fit, N_DRAWS)
  rnd    <- generate_common_randomness(length(params$sigma_mu), n_future)
  
  # Baseline
  base_logY <- simulate_intervention(params, rnd, 0, 0, start_idx = n_future + 1, n_future)
  base_2035 <- exp(base_logY[, n_future])
  counter <- counter + 1
  step6_rows[[counter]] <- tibble(Sector = pretty, Reduction_pct = 0, Gamma = 0, Damping = 0, Scenario = "Baseline",
                                  Y2035_Median = median(base_2035),Y2035_Q025 = quantile(base_2035, 0.025),
                                  Y2035_Q975 = quantile(base_2035, 0.975),
                                  P_exc = mean(base_logY[, n_future] > target_log))
  
  for (r in reduction_grid) {
    for (damp in damping_grid) {
      if (r == 0 && damp == 0) next
      if (r == 0 && damp > 0)  next
      scenario <- if (damp == 0) "Shock only" else "Shock + damping"
      logY <- simulate_intervention(params, rnd, r, damp, intervention_start_index, n_future)
      y2035 <- exp(logY[, n_future])
      counter <- counter + 1
      step6_rows[[counter]] <- tibble(Sector = pretty, Reduction_pct = 100 * r, Gamma = log(1 - r),
                                      Damping = damp, Scenario = scenario, Y2035_Median = median(y2035),
                                      Y2035_Q025 = quantile(y2035, 0.025),Y2035_Q975 = quantile(y2035, 0.975),
                                      P_exc = mean(logY[, n_future] > target_log))
    }
  }
}

step6_table <- bind_rows(step6_rows) %>%
  mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS), Damping = as.numeric(Damping)) %>%
  arrange(Sector, Reduction_pct, Damping)

shock_only_table <- step6_table %>%filter(Scenario == "Shock only") %>%
  select(Sector, Reduction_pct, Gamma, Y2035_Median, Y2035_Q025, Y2035_Q975, P_exc)

shock_damping_table <- step6_table %>% filter(Scenario == "Shock + damping") %>%
  select(Sector, Reduction_pct, Damping, Gamma, Y2035_Median,Y2035_Q025, Y2035_Q975, P_exc)

monotonicity_audit <- step6_table %>% filter(Scenario == "Shock only") %>% arrange(Sector, Reduction_pct) %>%
  group_by(Sector) %>% summarise(All_nonincreasing = all(diff(P_exc) <= 1e-8),
                                 Max_positive_change = max(diff(P_exc), na.rm = TRUE), .groups = "drop")

step6_results <- list(metadata = list(intervention_year = intervention_year,reduction_grid = reduction_grid,
                                      gamma_definition = "gamma = log(1 - reduction_fraction)",
                                      damping_grid = damping_grid,
                                      shock_definition = "Permanent level shift from 2025 onward",
                                      damping_definition = "Proportional reduction in post-intervention latent innovation scale",
                                      interpretation = paste("Conditional counterfactual dose-response.",
                                                             "Not a causally identified policy effect.",
                                                             "No threshold assumed a priori.")),
                      results= step6_table,shock_only= shock_only_table,shock_damping= shock_damping_table,
                      monotonicity= monotonicity_audit)
saveRDS(step6_results, file.path(OUTPUT_DIR, "Step6_intervention.rds"))
write_csv(step6_table,          file.path(OUTPUT_DIR, "Step6_Dose_Response.csv"))
write_csv(monotonicity_audit,   file.path(OUTPUT_DIR, "Step6_Monotonicity_Audit.csv"))

cat("06_intervention.R complete\n")
