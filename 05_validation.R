# =============================================================================
# 05_validation.R
# 05A Dynamic misspecification DGP calibration
# 05B Correlated-predictor DGP (optional)
# =============================================================================
message("=== 05_validation.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
# -----------------------------------------------------------------------------
# 5A Representative parameters (代表参数)
# -----------------------------------------------------------------------------
extract_representative_parameters <- function(fit) {
  post <- rstan::extract(fit, pars = c("s_mu", "s_Y", "mu_trend"))
  mu_trend <- post$mu_trend
  list(sigma_mu = median(post$s_mu), sigma_Y  = median(post$s_Y),
       mu_T1= median(mu_trend[, ncol(mu_trend) - 1]),mu_T= median(mu_trend[, ncol(mu_trend)]))
}

simulation_sectors <- c("coal", "oil")
representative_parameters <- purrr::map(simulation_sectors, function(s) {
  extract_representative_parameters(step2_results$fits[[s]][["Reduced_SORW"]])
})
names(representative_parameters) <- simulation_sectors
# -----------------------------------------------------------------------------
# 5B DGP generators
# -----------------------------------------------------------------------------
simulate_sorw_data <- function(T, mu_T1, mu_T, sigma_mu, sigma_Y) {
  mu <- numeric(T); mu[1] <- mu_T1; mu[2] <- mu_T
  if (T >= 3) for (t in 3:T) mu[t] <- 2*mu[t-1] - mu[t-2] + rnorm(1, 0, sigma_mu)
  y_log <- mu + rnorm(T, 0, sigma_Y)
  list(mu = mu, y_log = y_log, Y = exp(y_log))
}

simulate_forw_data <- function(T, mu_T, sigma_mu, sigma_Y) {
  mu <- numeric(T); mu[1] <- mu_T
  if (T >= 2) for (t in 2:T) mu[t] <- mu[t-1] + rnorm(1, 0, sigma_mu)
  y_log <- mu + rnorm(T, 0, sigma_Y)
  list(mu = mu, y_log = y_log, Y = exp(y_log))
}

simulate_ar1_data <- function(T, mu_start, rho = 0.7, sigma_mu, sigma_Y) {
  long_run_mean <- mu_start
  sigma_eta <- sigma_mu * sqrt(1 - rho^2)
  mu <- numeric(T); y_log <- numeric(T)
  mu[1] <- long_run_mean + rnorm(1, 0, sigma_mu)
  if (T >= 2) for (t in 2:T)
    mu[t] <- long_run_mean * (1 - rho) + rho * mu[t-1] + rnorm(1, 0, sigma_eta)
  y_log <- mu + rnorm(T, 0, sigma_Y)
  list(mu = mu, y_log = y_log, Y = exp(y_log),long_run_mean = long_run_mean, rho = rho, sigma_eta = sigma_eta)
}

simulate_structural_break_data <- function(T, mu_start, sigma_mu, sigma_Y,break_t = 8, break_size = NULL) {
  if (is.null(break_size)) break_size <- sigma_mu
  mu <- numeric(T); mu[1] <- mu_start
  if (T >= 2) for (t in 2:T) {
    mu[t] <- mu[t-1] + rnorm(1, 0, sigma_mu)
    if (t == break_t) mu[t] <- mu[t] + break_size
  }
  y_log <- mu + rnorm(T, 0, sigma_Y)
  list(mu = mu, y_log = y_log, Y = exp(y_log),break_t = break_t, break_size = break_size)
}
# -----------------------------------------------------------------------------
# 5C True DGP P_exc
# -----------------------------------------------------------------------------
true_p_exc_forw <- function(mu_T, sigma_mu, sigma_Y, Y_T,
                            target_ratio = CFG$target_ratio, horizon = 12) {
  threshold_log <- log(Y_T) + log(target_ratio)
  pred_sd <- sqrt(horizon * sigma_mu^2 + sigma_Y^2)
  1 - pnorm((threshold_log - mu_T) / pred_sd)
}

true_p_exc_sorw <- function(mu_T1, mu_T, sigma_mu, sigma_Y, Y_T,
                            target_ratio = CFG$target_ratio, horizon = 12) {
  threshold_log <- log(Y_T) + log(target_ratio)
  slope <- mu_T - mu_T1
  pred_mean <- mu_T + horizon * slope
  innov_var <- sigma_mu^2 * sum(seq_len(horizon)^2)
  pred_sd <- sqrt(innov_var + sigma_Y^2)
  1 - pnorm((threshold_log - pred_mean) / pred_sd)
}

true_p_exc_ar1 <- function(mu_T, long_run_mean, rho, sigma_eta, sigma_Y, Y_T,
                           target_ratio = CFG$target_ratio, horizon = 12) {
  threshold_log <- log(Y_T) + log(target_ratio)
  pred_mean <- long_run_mean + rho^horizon * (mu_T - long_run_mean)
  latent_var <- sigma_eta^2 * (1 - rho^(2*horizon)) / (1 - rho^2)
  pred_sd <- sqrt(latent_var + sigma_Y^2)
  1 - pnorm((threshold_log - pred_mean) / pred_sd)
}
# -----------------------------------------------------------------------------
# 5D Exact Grid Design Matrix & Posterior Solver
# -----------------------------------------------------------------------------
make_sorw_design <- function(T = 14) {
  Az <- matrix(0, T, T - 2)
  for (j in 1:(T - 2)) for (tau in (j + 2):T) Az[tau, j] <- tau - j - 1
  list(Az = Az, tv = 0:(T - 1))
}
SORW_DESIGN <- make_sorw_design(14)

fit_sorw_grid <- function(Y, covid = rep(0, length(Y)),shift = log((CFG$prior_alpha / CFG$prior_beta) / 0.5),
                          n_grid = 71, k_sd = 5, horizon = 12, target_ratio = CFG$target_ratio) {
  y <- log(Y); T <- length(y); Az <- SORW_DESIGN$Az; tv <- SORW_DESIGN$tv
  m_mu <- -2.0 + shift; sd_mu <- 0.7
  m_Y  <- -2.3 + shift; sd_Y  <- 0.6
  u_mu <- seq(m_mu - k_sd * sd_mu, m_mu + k_sd * sd_mu, length.out = n_grid)
  u_Y  <- seq(m_Y  - k_sd * sd_Y,  m_Y  + k_sd * sd_Y,  length.out = n_grid)
  base <- matrix(25, T, T) + tcrossprod(tv) + tcrossprod(covid)
  Q <- tcrossprod(Az); I_T <- diag(T)
  logw <- matrix(-Inf, n_grid, n_grid)
  
  for (i in seq_len(n_grid)) {
    S_i <- base + exp(2 * u_mu[i]) * Q
    for (k in seq_len(n_grid)) {
      R <- tryCatch(chol(S_i + exp(2 * u_Y[k]) * I_T), error = function(e) NULL)
      if (is.null(R)) next
      z  <- backsolve(R, y, transpose = TRUE)
      ll <- -0.5 * sum(z^2) - sum(log(diag(R))) - 0.5 * T * log(2 * pi)
      logw[i, k] <- ll + dnorm(u_mu[i], m_mu, sd_mu, log = TRUE) + dnorm(u_Y[k], m_Y, sd_Y, log = TRUE)
    }
  }
  w <- exp(logw - max(logw)); w <- w / sum(w)
  
  thr   <- y[T] + log(target_ratio)
  cells <- which(w > 1e-7, arr.ind = TRUE)
  pv0   <- c(25, 1, rep(1, T - 2), 1)
  h_var <- sum(seq_len(horizon)^2)
  p_cells <- numeric(nrow(cells))
  
  for (r in seq_len(nrow(cells))) {
    s1 <- exp(u_mu[cells[r, 1]]); s2 <- exp(u_Y[cells[r, 2]])
    A  <- cbind(1, tv, s1 * Az, covid)
    V  <- solve(diag(1 / pv0) + crossprod(A) / s2^2)
    m  <- V %*% crossprod(A, y) / s2^2
    aT  <- c(1, T - 1, s1 * Az[T, ],     0)
    aT1 <- c(1, T - 2, s1 * Az[T - 1, ], 0)
    b   <- aT + horizon * (aT - aT1)
    pv  <- drop(t(b) %*% V %*% b) + h_var * s1^2 + s2^2
    p_cells[r] <- 1 - pnorm((thr - sum(b * m)) / sqrt(pv))
  }
  wm <- rowSums(w); wy <- colSums(w)
  list(P_exc     = sum(w[cells] * p_cells) / sum(w[cells]),
       edge_mass = sum(w[c(1, n_grid), ]) + sum(w[, c(1, n_grid)]),
       s_mu_med  = exp(u_mu[which(cumsum(wm) >= 0.5)[1]]),
       s_Y_med   = exp(u_Y[which(cumsum(wy) >= 0.5)[1]]))
}

# -----------------------------------------------------------------------------
# 5E Main simulation loop
# -----------------------------------------------------------------------------
N_SIM <- 200
SIM_SEED <- 20260909
simulation_dgps <- c("SORW_correct", "FORW_misspecified","AR1_misspecified", "StructuralBreak_misspecified")
simulation_results <- list()
counter <- 0
total_fits <- length(simulation_sectors) * length(simulation_dgps) * N_SIM
pb <- utils::txtProgressBar(min = 0, max = total_fits, style = 3)
start_time <- Sys.time()
for (sector in simulation_sectors) {
  p <- representative_parameters[[sector]]
  for (dgp in simulation_dgps) {
    cat("\n\n--- Starting:", SECTOR_LABELS[[sector]], "| DGP:", dgp, "---\n")
    for (rep in seq_len(N_SIM)) {
      counter <- counter + 1
      seed_i <- SIM_SEED + counter
      set.seed(seed_i)
      # --- Simulate synthetic dataset ---
      sim_data <- switch(dgp,SORW_correct = simulate_sorw_data(14, p$mu_T1, p$mu_T, p$sigma_mu, p$sigma_Y),
                         FORW_misspecified = simulate_forw_data(14, p$mu_T, p$sigma_mu, p$sigma_Y),
                         AR1_misspecified = simulate_ar1_data(14, p$mu_T, 0.7, p$sigma_mu, p$sigma_Y),
                         StructuralBreak_misspecified =simulate_structural_break_data(14, p$mu_T, p$sigma_mu, p$sigma_Y))
      # --- True DGP P_exc ---
      true_p <- switch(dgp,SORW_correct = true_p_exc_sorw(sim_data$mu[length(sim_data$mu) - 1],
                                                          tail(sim_data$mu, 1), p$sigma_mu, p$sigma_Y,tail(sim_data$Y, 1)),
                       FORW_misspecified = true_p_exc_forw(tail(sim_data$mu, 1), p$sigma_mu, p$sigma_Y,tail(sim_data$Y, 1)),
                       AR1_misspecified = true_p_exc_ar1(tail(sim_data$mu, 1), sim_data$long_run_mean,0.7, sim_data$sigma_eta, 
                                                         p$sigma_Y,tail(sim_data$Y, 1)),
                       StructuralBreak_misspecified =true_p_exc_forw(tail(sim_data$mu, 1), p$sigma_mu, p$sigma_Y,
                                                                     tail(sim_data$Y, 1)))
      # --- Exact posterior via grid  ---
      g <- fit_sorw_grid(sim_data$Y)
      simulation_results[[counter]] <- tibble(Sector = SECTOR_LABELS[[sector]], DGP = dgp, Replicate = rep,
                                              Fit_OK = g$edge_mass < 0.01,
                                              True_P_exc = true_p, Fitted_P_exc = g$P_exc,
                                              P_exc_Error = g$P_exc - true_p, Max_Rhat = 1.00, Divergences = 0)
      utils::setTxtProgressBar(pb, counter)
      # --- Progress message every 25 fits ---
      if (counter %% 25 == 0 || counter == 1 || counter == total_fits) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        rate <- counter / elapsed
        remaining <- total_fits - counter
        eta_sec <- ifelse(rate > 0, remaining / rate, NA_real_)
        
        cat("\n[Progress] ", counter, "/", total_fits,
            " (", round(100 * counter / total_fits, 1), "%)",
            " | elapsed = ", round(elapsed / 60, 1), " min",
            " | rate = ", round(rate, 2), " fits/min",
            " | ETA ≈ ", round(eta_sec / 60, 1), " min\n", sep = "")
      }
      
      # --- Checkpoint every 100 fits ---
      # if (counter %% 100 == 0) {
      #   simulation_checkpoint <- bind_rows(simulation_results)
      #   saveRDS(simulation_checkpoint,file.path(OUTPUT_DIR,paste0("Step5_simulation_checkpoint_", counter, ".rds")))
      #   write_csv(simulation_checkpoint,file.path(OUTPUT_DIR,paste0("Step5_simulation_checkpoint_", counter, ".csv")))
      #   cat("[Checkpoint saved at ", counter, " fits]\n", sep = "")
      # }
    }
  }
}

close(pb)
total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n\nSimulation complete.", "\nTotal fits: ", counter,
    "\nTotal elapsed time: ", round(total_elapsed, 1), " min\n", sep = "")
# -----------------------------------------------------------------------------
# Combine results
# -----------------------------------------------------------------------------
simulation_df <- bind_rows(simulation_results)
simulation_df %>% group_by(Sector, DGP) %>% summarise(fail_rate = mean(!Fit_OK),
                                                      true_pexc_failed  = mean(True_P_exc[!Fit_OK]),
                                                      true_pexc_passed  = mean(True_P_exc[Fit_OK]),.groups = "drop")
# -----------------------------------------------------------------------------
# Calibration summary
# -----------------------------------------------------------------------------
calibration_summary <- simulation_df %>%filter(Fit_OK) %>%group_by(Sector, DGP) %>%
  summarise(N = n(),True_P_exc_Mean = mean(True_P_exc, na.rm = TRUE),
            Fitted_P_exc_Mean = mean(Fitted_P_exc, na.rm = TRUE),Mean_Error = mean(P_exc_Error, na.rm = TRUE),
            RMSE = sqrt(mean(P_exc_Error^2, na.rm = TRUE)),.groups = "drop")
# -----------------------------------------------------------------------------
# Save final results
# -----------------------------------------------------------------------------
step5_validation <- list(
  metadata = list(N_SIM = N_SIM, DGPs = simulation_dgps, sectors = simulation_sectors,
                  total_fits = total_fits, fitting_model = "Reduced SORW",
                  interpretation = "Calibration/misspecification sensitivity, not robustness proof."),
  representative_parameters = representative_parameters,raw = simulation_df, summary = calibration_summary)

saveRDS(step5_validation, file.path(OUTPUT_DIR, "Step5_validation.rds"))
write_csv(simulation_df,       file.path(OUTPUT_DIR, "Step5_Simulation_Raw.csv"))
write_csv(calibration_summary, file.path(OUTPUT_DIR, "Step5_Simulation_Summary.csv"))

cat("\n05_validation.R complete\n")