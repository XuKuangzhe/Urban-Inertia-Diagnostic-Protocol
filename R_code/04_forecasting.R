# =============================================================================
# 04_forecasting.R
# PRIMARY P_exc is produced solely by analytical_p_exc().
# (PRIMARY P_exc 由 analytical_p_exc() 唯一产生)
# Monte Carlo is used only for QC and fan charts.
# (Monte Carlo 仅用于 QC 和 fan chart)
# Input (输入)：Step2_primary_results.rds
# Output (输出)：Step4_forecast.rds
# =============================================================================
message("=== 04_forecasting.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
# -----------------------------------------------------------------------------
# 4A Analytical (解析) P_exc
# -----------------------------------------------------------------------------
analytic_p_exc <- function(fit, Y_T, model = c("SORW", "FORW", "DRIFT"),target_ratio = CFG$target_ratio, horizons = 1:12) {
  model <- match.arg(model)
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y", if (model == "DRIFT") "drift"))
  mu_trend <- post$mu_trend
  mu_T1 <- as.numeric(mu_trend[, ncol(mu_trend) - 1])
  mu_T  <- as.numeric(mu_trend[, ncol(mu_trend)])
  sigma_mu <- as.numeric(post$s_mu); sigma_Y <- as.numeric(post$s_Y)
  threshold_log <- log(Y_T) + log(target_ratio)
  
  purrr::map_dfr(horizons, function(h) {
    if (model == "FORW") {
      pred_mean <- mu_T
      pred_var  <- h * sigma_mu^2 + sigma_Y^2
    } else if (model == "DRIFT") {
      pred_mean <- mu_T + h * as.numeric(post$drift)
      pred_var  <- h * sigma_mu^2 + sigma_Y^2
    } else {
      slope     <- mu_T - mu_T1
      pred_mean <- mu_T + h * slope
      innov_var <- sum(seq_len(h)^2) * sigma_mu^2
      pred_var  <- innov_var + sigma_Y^2
    }
    cond_p <- 1 - pnorm(threshold_log, mean = pred_mean, sd = sqrt(pred_var))
    tibble(Horizon = h, Target_Ratio = target_ratio, P_exc = mean(cond_p))
  })
}

analytic_p_exc_table <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  pretty <- SECTOR_LABELS[[sector]]
  Y_T <- tail(as.numeric(D[[sector]]), 1)
  fit_sorw  <- step2_results$fits[[sector]]$Reduced_SORW
  fit_forw  <- step2_results$fits[[sector]]$Reduced_FORW
  fit_drift <- step2_results$fits[[sector]]$Reduced_DriftFORW
  
  bind_rows(analytic_p_exc(fit_sorw,  Y_T, "SORW",  CFG$target_ratio, 1:12) %>% mutate(Sector = pretty, Model = "SORW"),
            analytic_p_exc(fit_forw,  Y_T, "FORW",  CFG$target_ratio, 1:12) %>% mutate(Sector = pretty, Model = "FORW"),
            analytic_p_exc(fit_drift, Y_T, "DRIFT", CFG$target_ratio, 1:12) %>% mutate(Sector = pretty, Model = "DRIFT"))
})

analytic_2035_p_exc <- analytic_p_exc_table %>% filter(Horizon == 12) %>% select(Sector, Model, P_exc) %>%
  pivot_wider(names_from = Model, values_from = P_exc, names_prefix = "P_exc_") %>%
  mutate(Diff_SORW_minus_FORW  = P_exc_SORW - P_exc_FORW, Diff_SORW_minus_DRIFT = P_exc_SORW - P_exc_DRIFT)

# -----------------------------------------------------------------------------
# 4B Horizon dependence
# -----------------------------------------------------------------------------
# Already covered in analytic_p_exc_table; used directly for plotting.
# (analytic_p_exc_table 已覆盖，直接用于绘图)

# -----------------------------------------------------------------------------
# 4C Target-ratio sensitivity
# -----------------------------------------------------------------------------
target_ratio_grid <- seq(0.50, 1.20, by = 0.02)
target_ratio_results <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  pretty <- SECTOR_LABELS[[sector]]
  Y_T <- tail(as.numeric(D[[sector]]), 1)
  fit_sorw <- step2_results$fits[[sector]]$Reduced_SORW
  fit_forw <- step2_results$fits[[sector]]$Reduced_FORW
  
  bind_rows(
    purrr::map_dfr(target_ratio_grid, function(r) {
      analytic_p_exc(fit_sorw, Y_T, "SORW", r, 12) %>%
        mutate(Sector = pretty, Model = "SORW")
    }),
    purrr::map_dfr(target_ratio_grid, function(r) {
      analytic_p_exc(fit_forw, Y_T, "FORW", r, 12) %>%
        mutate(Sector = pretty, Model = "FORW")
    })
  )
})

# -----------------------------------------------------------------------------
# 4D Monte Carlo forecast (仅 QC + fan chart)
# -----------------------------------------------------------------------------
simulate_sorw_forecast <- function(fit, n_future = length(CFG$future_years)) {
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y"))
  n_draws <- length(post$s_mu)
  mu_draws <- post$mu_trend
  mu_T1 <- mu_draws[, ncol(mu_draws) - 1]
  mu_T  <- mu_draws[, ncol(mu_draws)]
  sigma_mu <- as.numeric(post$s_mu)
  sigma_Y  <- as.numeric(post$s_Y)
  
  future_logY <- matrix(NA_real_, n_draws, n_future)
  prev <- mu_T1; curr <- mu_T
  for (h in seq_len(n_future)) {
    next_mu   <- 2 * curr - prev + sigma_mu * rnorm(n_draws)
    next_logY <- next_mu + sigma_Y * rnorm(n_draws)
    future_logY[, h] <- next_logY
    prev <- curr; curr <- next_mu
  }
  future_logY
}

simulate_forw_forecast <- function(fit, n_future = length(CFG$future_years)) {
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y"))
  n_draws <- length(post$s_mu)
  mu_draws <- post$mu_trend
  curr <- mu_draws[, ncol(mu_draws)]
  sigma_mu <- as.numeric(post$s_mu)
  sigma_Y  <- as.numeric(post$s_Y)
  
  future_logY <- matrix(NA_real_, n_draws, n_future)
  for (h in seq_len(n_future)) {
    next_mu   <- curr + sigma_mu * rnorm(n_draws)
    next_logY <- next_mu + sigma_Y * rnorm(n_draws)
    future_logY[, h] <- next_logY
    curr <- next_mu
  }
  future_logY
}

summarize_forecast <- function(logY_draws, years) {
  Y_draws <- exp(logY_draws)
  tibble(Year = years, Mean = apply(Y_draws, 2, mean), Median = apply(Y_draws, 2, median),
         Lower_95 = apply(Y_draws, 2, quantile, probs = 0.025),
         Upper_95 = apply(Y_draws, 2, quantile, probs = 0.975),
         Log_Mean = apply(logY_draws, 2, mean), Log_Median = apply(logY_draws, 2, median),
         Log_Lower_95 = apply(logY_draws, 2, quantile, probs = 0.025),
         Log_Upper_95 = apply(logY_draws, 2, quantile, probs = 0.975))
}

set.seed(CFG$seed)
forecast_long <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  pretty <- SECTOR_LABELS[[sector]]
  fit_sorw <- step2_results$fits[[sector]]$Reduced_SORW
  fit_forw <- step2_results$fits[[sector]]$Reduced_FORW
  
  logY_sorw <- simulate_sorw_forecast(fit_sorw)
  logY_forw <- simulate_forw_forecast(fit_forw)
  
  bind_rows(summarize_forecast(logY_sorw, CFG$future_years) %>% mutate(Sector = pretty, Model = "SORW"),
            summarize_forecast(logY_forw, CFG$future_years) %>% mutate(Sector = pretty, Model = "FORW"))
})

# -----------------------------------------------------------------------------
# 4E MC vs analytical QC
# -----------------------------------------------------------------------------
mc_p_exc <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  pretty <- SECTOR_LABELS[[sector]]
  Y_2023 <- tail(as.numeric(D[[sector]]), 1)
  target_log <- log(Y_2023) + log(CFG$target_ratio)
  fit_sorw <- step2_results$fits[[sector]]$Reduced_SORW
  fit_forw <- step2_results$fits[[sector]]$Reduced_FORW
  
  set.seed(CFG$seed + 1)
  logY_sorw <- simulate_sorw_forecast(fit_sorw)
  set.seed(CFG$seed + 1)
  logY_forw <- simulate_forw_forecast(fit_forw)
  
  tibble(Sector = pretty,
         P_exc_SORW_MC = mean(logY_sorw[, ncol(logY_sorw)] > target_log),
         P_exc_FORW_MC = mean(logY_forw[, ncol(logY_forw)] > target_log))
})

mc_vs_analytic <- mc_p_exc %>% left_join(analytic_2035_p_exc %>%select(Sector, P_exc_SORW_analytic = P_exc_SORW,
                                                                       P_exc_FORW_analytic = P_exc_FORW),
                                         by = "Sector") %>%mutate(SORW_diff = P_exc_SORW_MC - P_exc_SORW_analytic,
                                                                  FORW_diff = P_exc_FORW_MC - P_exc_FORW_analytic)

# =============================================================================
# save (保存)
# =============================================================================
step4_forecast <- list(
  metadata = list(forecast_years = CFG$future_years,
                  target_definition = sprintf("Y2035 <= %.2f * Y2023", CFG$target_ratio),
                  primary_P_exc_source = "analytic_p_exc()", MC_role = "QC and fan chart only"),
  analytic_p_exc_table    = analytic_p_exc_table,
  analytic_2035_p_exc     = analytic_2035_p_exc,
  target_ratio_results    = target_ratio_results,
  forecast_long           = forecast_long,
  mc_vs_analytic          = mc_vs_analytic
)
saveRDS(step4_forecast, file.path(OUTPUT_DIR, "Step4_forecast.rds"))
write_csv(analytic_p_exc_table,  file.path(OUTPUT_DIR, "Step4_Analytic_Pexc_Horizon.csv"))
write_csv(analytic_2035_p_exc,   file.path(OUTPUT_DIR, "Step4_Analytic_Pexc_2035.csv"))
write_csv(target_ratio_results,  file.path(OUTPUT_DIR, "Step4_Pexc_TargetRatio.csv"))
write_csv(mc_vs_analytic,        file.path(OUTPUT_DIR, "Step4_MC_vs_Analytic_QC.csv"))

cat("04_forecasting.R complete\n")
