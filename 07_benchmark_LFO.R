# =============================================================================
# 07_benchmark_LFO.R
# 07A 4-way dynamic benchmark: SORW / FORW / DriftFORW / LLT
# 07B LFO validation
# =============================================================================
message("=== 07_benchmark_LFO.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
MODELS$Reduced_DriftFORW <- stan_model(file.path(STAN_DIR, "SSM_Reduced_DriftFORW_v1.stan"))
MODELS$Reduced_LLT       <- stan_model(file.path(STAN_DIR, "SSM_Reduced_LLT_v1.stan"))
# -----------------------------------------------------------------------------
# 7A 4-way comparison
# -----------------------------------------------------------------------------
fit_benchmark_model <- function(model, data, seed, model_name, sector_name,iter = CFG$iter,
                                adapt_delta_grid = c(CFG$adapt_delta, 0.995, 0.999),max_treedepth = CFG$max_treedepth) {
  adapt_delta_grid <- unique(adapt_delta_grid)
  for (attempt in seq_along(adapt_delta_grid)) {
    delta_i <- adapt_delta_grid[attempt]
    cat("\n  -> ", sector_name, " | ", model_name, " | attempt ", attempt, "/",
        length(adapt_delta_grid), " | adapt_delta = ", delta_i, "\n", sep = "")
    fit_start <- Sys.time()
    fit <- tryCatch(sampling(model, data = data, iter = iter, seed = seed,
                             control = list(adapt_delta = delta_i, max_treedepth = max_treedepth),refresh = 0),
                    error = function(e) {
                      cat("     Sampling error: ", conditionMessage(e), "\n", sep = "")
                      NULL
                    }
    )
    fit_elapsed <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
    if (is.null(fit)) next
    diag <- diagnose_stan_fit(fit, model_name = model_name, sector = sector_name)
    all_ok <- diag$Rhat_OK && diag$Divergence_OK && diag$Treedepth_OK && diag$BFMI_OK
    cat("     Finished in ", round(fit_elapsed, 2), " min | divergences = ",
        diag$Divergences, " | max Rhat = ", round(diag$Max_Rhat, 5),
        " | BFMI = ", round(diag$Min_BFMI, 3), "\n", sep = "")
    if (all_ok) {
      return(list(fit = fit, diagnostics = diag, metadata = list(seed = seed, iter = iter, adapt_delta_used = delta_i,
                                                                 retry_number = attempt - 1, elapsed_minutes = fit_elapsed)))
    }
    cat("     Diagnostic failure -> retrying if possible.\n")
  }
  stop(paste0("Benchmark model failed diagnostics after all retries: ", sector_name, " / ", model_name))
}

fit_dynamic_benchmarks <- function(Y, seed, sector_name) {
  data_i <- make_reduced_data(Y)
  fit_drift <- fit_benchmark_model(model = MODELS$Reduced_DriftFORW, data = data_i,seed = seed + 1, model_name = "Drift-FORW",
                                   sector_name = sector_name)
  fit_llt <- fit_benchmark_model(model = MODELS$Reduced_LLT, data = data_i,seed = seed + 2, model_name = "LLT",
                                 sector_name = sector_name)
  list(DriftFORW = fit_drift, LLT = fit_llt)
}

dynamic_comparison_results <- list()
n_sectors <- length(SECTOR_CODES)
total_new_fits <- n_sectors * 2
progress_counter <- 0
start_time <- Sys.time()
pb <- utils::txtProgressBar(min = 0, max = total_new_fits, style = 3)
for (sector_index in seq_along(SECTOR_CODES)) {
  sector <- SECTOR_CODES[sector_index]
  pretty <- SECTOR_LABELS[[sector]]
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Sector ", sector_index, "/", n_sectors, ": ", pretty, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  Y <- as.numeric(D[[sector]])
  sector_seed <- 9000 + 100 * sector_index
  # --- Fit Drift-FORW ---
  drift_result <- fit_benchmark_model(model = MODELS$Reduced_DriftFORW,data = make_reduced_data(Y),
                                      seed = sector_seed + 1,model_name = "Drift-FORW",sector_name = pretty)
  progress_counter <- progress_counter + 1
  utils::setTxtProgressBar(pb, progress_counter)
  # --- Fit LLT ---
  llt_result <- fit_benchmark_model(model = MODELS$Reduced_LLT,data = make_reduced_data(Y),seed = sector_seed + 2,
                                    model_name = "LLT",sector_name = pretty)
  progress_counter <- progress_counter + 1
  utils::setTxtProgressBar(pb, progress_counter)
  # --- Construct LOO objects ---
  loo_objects_dyn <- list(SORW = step2_results$loo[[sector]]$Reduced_SORW,
                          FORW = step2_results$loo[[sector]]$Reduced_FORW,
                          DriftFORW = loo(extract_log_lik(drift_result$fit, parameter_name = "log_lik",merge_chains = FALSE)),
                          LLT = loo(extract_log_lik(llt_result$fit, parameter_name = "log_lik", merge_chains = FALSE)))
  # --- Four-way comparison ---
  comparison <- loo_compare(loo_objects_dyn)
  dynamic_comparison_results[[sector]] <- list(fits = list(DriftFORW = drift_result$fit, LLT = llt_result$fit),
                                               diagnostics = list(DriftFORW = drift_result$diagnostics,
                                                                  LLT = llt_result$diagnostics),
                                               fit_metadata = list(DriftFORW = drift_result$metadata,
                                                                   LLT = llt_result$metadata),
                                               loo = loo_objects_dyn,comparison = comparison)
  # --- Print sector result ---
  cat("\n====================================================\n", pretty,
      "\n====================================================\n")
  print(comparison)
  # --- Progress / elapsed time ---
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  fits_done <- progress_counter
  fit_rate <- ifelse(elapsed > 0, fits_done / elapsed, NA_real_)
  remaining <- total_new_fits - fits_done
  eta_minutes <- ifelse(is.finite(fit_rate) && fit_rate > 0, remaining / fit_rate, NA_real_)
  
  cat("\n[07A Progress] ", fits_done, "/", total_new_fits, " fits (",
      round(100 * fits_done / total_new_fits, 1), "%)",
      " | elapsed = ", round(elapsed, 1), " min",
      " | ETA ≈ ", round(eta_minutes, 1), " min\n", sep = "")
  # --- Checkpoint after each sector ---
  # saveRDS(dynamic_comparison_results,file.path(OUTPUT_DIR, "Step7_dynamic_benchmark_checkpoint.rds"))
  # cat("Checkpoint saved after sector: ", pretty, "\n", sep = "")
}

close(pb)

step7_benchmark <- list(metadata = list(models = c("Reduced SORW", "Reduced FORW", "Drift-FORW", "LLT"),
                                        sectors = SECTOR_CODES,total_new_fits = total_new_fits,
                                        interpretation ="Four-way predictive benchmark comparison; not a causal model ranking."),
                        results = dynamic_comparison_results)
saveRDS(step7_benchmark,file.path(OUTPUT_DIR, "Step7_dynamic_benchmark_results.rds"))

# -----------------------------------------------------------------------------
# 7B Predictive density functions
# -----------------------------------------------------------------------------
predictive_log_density_forw <- function(fit, y_log, horizon, covid_target) {
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y", "beta_covid"))
  mu_T <- post$mu_trend[, ncol(post$mu_trend)]
  mean_h <- mu_T + post$beta_covid * covid_target
  var_h  <- horizon * post$s_mu^2 + post$s_Y^2
  ld <- dnorm(y_log, mean = mean_h, sd = sqrt(var_h), log = TRUE)
  log_mean_exp(ld)
}

predictive_log_density_drift <- function(fit, y_log, horizon, covid_target) {
  post <- rstan::extract(fit, pars = c("mu_trend", "drift", "s_mu", "s_Y", "beta_covid"))
  mu_T <- post$mu_trend[, ncol(post$mu_trend)]
  mean_h <- mu_T + horizon * as.numeric(post$drift) + post$beta_covid * covid_target
  var_h  <- horizon * post$s_mu^2 + post$s_Y^2
  ld <- dnorm(y_log, mean = mean_h, sd = sqrt(var_h), log = TRUE)
  log_mean_exp(ld)
}

predictive_log_density_sorw <- function(fit, y_log, horizon, covid_target) {
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y", "beta_covid"))
  mu_T1 <- post$mu_trend[, ncol(post$mu_trend) - 1]
  mu_T  <- post$mu_trend[, ncol(post$mu_trend)]
  slope <- mu_T - mu_T1
  mean_h <- mu_T + horizon * slope + post$beta_covid * covid_target
  var_h  <- sum(seq_len(horizon)^2) * post$s_mu^2 + post$s_Y^2
  ld <- dnorm(y_log, mean = mean_h, sd = sqrt(var_h), log = TRUE)
  log_mean_exp(ld)
}

# LLT: level innovation weight sum_{j=1}^h j^2;
#      slope innovation weight sum_{j=1}^{h-1} j^2
# (LLT：level innovation 权重 sum_{j=1}^h j^2；slope innovation 权重 sum_{j=1}^{h-1} j^2)
predictive_log_density_llt <- function(fit, y_log, horizon, covid_target) {
  post <- rstan::extract(fit, pars = c("level", "slope", "s_level", "s_slope", "s_Y", "beta_covid"))
  level_T <- post$level[, ncol(post$level)]
  slope_T <- post$slope[, ncol(post$slope)]
  mean_h  <- level_T + horizon * slope_T + post$beta_covid * covid_target
  w_level <- sum(seq_len(horizon)^2)
  w_slope <- if (horizon <= 1) 0 else sum(seq_len(horizon - 1)^2)
  var_h   <- w_level * post$s_level^2 + w_slope * post$s_slope^2 + post$s_Y^2
  ld <- dnorm(y_log, mean = mean_h, sd = sqrt(var_h), log = TRUE)
  log_mean_exp(ld)
}

# -----------------------------------------------------------------------------
# 7C LFO
# -----------------------------------------------------------------------------
lfo_origins  <- 2018:2021
lfo_horizons <- 1:2
fit_prefix_model <- function(model, Y, covid_dummy, seed, iter = 2000) {
  n <- length(Y)
  data_i <- list(T = n, Y = as.numeric(Y), covid_dummy = as.numeric(covid_dummy[seq_len(n)]),
                 prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta)
  sampling(model, data = data_i, iter = iter, seed = seed,
           control = list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth), refresh = 0)
}

tasks_per_sector <- sum(sapply(lfo_origins, function(origin) {
  sum((origin + lfo_horizons) %in% CFG$years)
}))
total_tasks <- length(SECTOR_CODES) * tasks_per_sector
lfo_results <- list()
counter <- 0

cat(sprintf("\n>>> Starting 7C LFO: %d sectors × %d evals = %d total tasks\n", 
            length(SECTOR_CODES), tasks_per_sector, total_tasks))
pb <- utils::txtProgressBar(min = 0, max = total_tasks, style = 3)
for (sector in SECTOR_CODES) {
  pretty <- SECTOR_LABELS[[sector]]
  Y_full <- as.numeric(D[[sector]])
  
  for (origin_year in lfo_origins) {
    train_idx <- which(CFG$years <= origin_year)
    for (h in lfo_horizons) {
      target_year <- origin_year + h
      if (!(target_year %in% CFG$years)) next
      
      counter <- counter + 1
      Y_train <- Y_full[train_idx]
      covid_target <- as.numeric(target_year %in% CFG$covid_years)
      
      fit_sorw  <- fit_prefix_model(MODELS$Reduced_SORW,      Y_train, COVID_PRIMARY, 11000 + counter)
      fit_forw  <- fit_prefix_model(MODELS$Reduced_FORW,      Y_train, COVID_PRIMARY, 12000 + counter)
      fit_drift <- fit_prefix_model(MODELS$Reduced_DriftFORW, Y_train, COVID_PRIMARY, 13000 + counter)
      
      fit_llt_data <- list(T = length(Y_train), Y = Y_train, covid_dummy = as.numeric(COVID_PRIMARY[seq_along(Y_train)]),
                           prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta)
      fit_llt <- sampling(MODELS$Reduced_LLT, data = fit_llt_data, iter = 2000, seed = 14000 + counter,
                          control = list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth), refresh = 0)
      
      y_log_target <- log(Y_full[which(CFG$years == target_year)])
      
      lfo_results[[counter]] <- tibble(
        Sector = pretty, Origin = origin_year, TargetYear = target_year, Horizon = h,
        COVID_Target = target_year %in% CFG$covid_years,
        SORW      = predictive_log_density_sorw(fit_sorw, y_log_target, h, covid_target),
        FORW      = predictive_log_density_forw(fit_forw, y_log_target, h, covid_target),
        DriftFORW = predictive_log_density_drift(fit_drift, y_log_target, h, covid_target),
        LLT       = predictive_log_density_llt(fit_llt, y_log_target, h, covid_target)
      )
      utils::setTxtProgressBar(pb, counter)
    }
  }
}
close(pb)

lfo_table <- bind_rows(lfo_results)
lfo_summary <- lfo_table %>% pivot_longer(cols = c(SORW, FORW, DriftFORW, LLT),
                                          names_to = "Model", values_to = "LogPredictiveDensity") %>%
  group_by(Sector, Model) %>% summarise(Mean_LPD = mean(LogPredictiveDensity),
                                        SD_LPD = sd(LogPredictiveDensity), N = n(), .groups = "drop")

step7_benchmark <- list(results = dynamic_comparison_results, lfo_table = lfo_table, lfo_summary = lfo_summary)
saveRDS(step7_benchmark, file.path(OUTPUT_DIR, "Step7_benchmark.rds"))
write_csv(lfo_table,   file.path(OUTPUT_DIR, "Step7_LFO_PredictiveScores.csv"))
write_csv(lfo_summary, file.path(OUTPUT_DIR, "Step7_LFO_Summary.csv"))

cat("07_benchmark_LFO.R complete\n")