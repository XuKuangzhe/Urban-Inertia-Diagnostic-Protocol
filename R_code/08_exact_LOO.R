# =============================================================================
# 08_exact_LOO.R
# Exact leave-one-out audit for k > 0.7 points
# Simplified triage: refit all k > 0.7; flag k > 1.0 as highest priority.
# (简化 triage：k > 0.7 全部 refit；k > 1.0 标 highest priority)
# input (输入)：Step2_primary_results.rds
# output (输出)：Step8_exact_LOO.rds
# =============================================================================
message("=== 08_exact_LOO.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
EXACT_ITER        <- 4000
EXACT_WARMUP      <- 2000
EXACT_CHAINS      <- 4
EXACT_SEED_BASE   <- 70000
# -----------------------------------------------------------------------------
# 8A Build v3 data (with use_obs)
# 8A 构造 v3 data（含 use_obs）
# -----------------------------------------------------------------------------
make_full_data_v3 <- function(Y, use_obs = rep(1L, length(CFG$years))) {
  list(T = length(Y), K = ncol(X_FULL), X = X_FULL, Y = as.numeric(Y), 
       covid_dummy = as.numeric(COVID_PRIMARY), use_obs = as.integer(use_obs),
       prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta,
       lasso_alpha = CFG$lasso_alpha, lasso_beta = CFG$lasso_beta)
}

make_reduced_data_v3 <- function(Y, use_obs = rep(1L, length(CFG$years))) {
  list(T = length(Y), Y = as.numeric(Y),covid_dummy = as.numeric(COVID_PRIMARY), use_obs = as.integer(use_obs),
       prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta)
}

step2_data_v3 <- list()
for (sector in SECTOR_CODES) {
  Y <- as.numeric(D[[sector]])
  step2_data_v3[[sector]] <- list(Full_SORW= make_full_data_v3(Y),Reduced_SORW= make_reduced_data_v3(Y),
                                  Reduced_FORW= make_reduced_data_v3(Y),Reduced_DriftFORW = make_reduced_data_v3(Y))
}

# -----------------------------------------------------------------------------
# 8B Exact LOO on one point
# -----------------------------------------------------------------------------
exact_loo_one_point <- function(stan_model_obj, base_data, held_out, year,seed = EXACT_SEED_BASE) {
  T <- base_data$T
  stopifnot(held_out >= 1, held_out <= T)
  loo_data <- base_data
  loo_data$use_obs <- rep(1L, T)
  loo_data$use_obs[held_out] <- 0L
  
  fit <- sampling(stan_model_obj, data = loo_data, chains = EXACT_CHAINS,iter = EXACT_ITER, warmup = EXACT_WARMUP, 
                  cores = EXACT_CHAINS,seed = seed + held_out,
                  control = list(adapt_delta = CFG$adapt_delta,max_treedepth = CFG$max_treedepth),refresh = 0)
  
  ll <- as.vector(rstan::extract(fit, pars = "log_lik", permuted = FALSE)[, , held_out])
  stopifnot(all(is.finite(ll)))
  lpd <- log_mean_exp(ll)
  
  summ <- summary(fit)$summary
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(held_out = held_out, year = year, lpd_exact = lpd,max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
       min_n_eff = min(summ[, "n_eff"], na.rm = TRUE),
       divergences = sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1))),
       treedepth_hits = sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= CFG$max_treedepth), numeric(1))))
}

# -----------------------------------------------------------------------------
# 8C Target (目标点)
# -----------------------------------------------------------------------------
high_k_targets <- step2_results$pareto_audit %>% filter(Pareto_k > 0.7)
n_targets <- nrow(high_k_targets)
exact_results <- list()

cat(sprintf("Exact LOO targets (k > 0.7): %d points\n", n_targets))
if (n_targets > 0) {
  pb <- utils::txtProgressBar(min = 0, max = n_targets, style = 3)
  
  for (i in seq_len(n_targets)) {
    row <- high_k_targets[i, ]
    sector <- row$Sector; model_name <- row$Model
    idx <- row$Index; yr <- row$Year
    v3_model <- MODELS[[paste0(model_name, "_v3")]]
    
    if (is.null(v3_model)) {
      warning(sprintf("No v3 model for %s; skipping %s %d", model_name, sector, yr))
      utils::setTxtProgressBar(pb, i)
      next
    }
    
    r <- exact_loo_one_point(
      v3_model, step2_data_v3[[sector]][[model_name]], held_out = idx, year = yr,
      seed = EXACT_SEED_BASE + 1000 * match(sector, SECTOR_CODES) + idx
    )
    
    exact_results[[i]] <- tibble(
      Sector = sector, Model = model_name, Year = yr, Index = idx, Exact_ELPD = r$lpd_exact,
      Max_Rhat = r$max_rhat, Min_n_eff = r$min_n_eff,
      Divergences = r$divergences, Treedepth_hits = r$treedepth_hits
    )
    
    utils::setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("\n>>> Exact LOO refitting complete.\n")
} else {
  cat(">>> No observations with Pareto-k > 0.7. Skipping exact LOO refits.\n")
}

exact_audit_results <- bind_rows(exact_results) %>%
  mutate(Div_Rate = Divergences / (EXACT_CHAINS * (EXACT_ITER - EXACT_WARMUP)),
         Convergence_OK = Max_Rhat < 1.01 & Min_n_eff > 400 & Div_Rate < 0.01)

if (nrow(exact_audit_results) > 0 && !all(exact_audit_results$Convergence_OK)) {
  warning("Some exact-LOO refits failed convergence checks.")
}

# -----------------------------------------------------------------------------
# 8D Hybrid comparison
# -----------------------------------------------------------------------------
build_hybrid_pointwise <- function(loo_obj, exact_rows) {
  pw <- as.data.frame(loo_obj$pointwise)
  if (nrow(exact_rows) > 0) {
    for (i in seq_len(nrow(exact_rows))) {
      pw$elpd_loo[exact_rows$Index[i]] <- exact_rows$Exact_ELPD[i]
    }
  }
  pw
}

compare_hybrid <- function(sector, model_1, model_2) {
  loo_1 <- step2_results$loo[[sector]][[model_1]]
  loo_2 <- step2_results$loo[[sector]][[model_2]]
  ex_1 <- exact_audit_results %>% filter(Sector == sector, Model == model_1)
  ex_2 <- exact_audit_results %>% filter(Sector == sector, Model == model_2)
  pw_1 <- build_hybrid_pointwise(loo_1, ex_1)
  pw_2 <- build_hybrid_pointwise(loo_2, ex_2)
  delta_i <- pw_1$elpd_loo - pw_2$elpd_loo
  tibble(Sector = SECTOR_LABELS[[sector]], Model_1 = model_1, Model_2 = model_2,ELPD_1 = sum(pw_1$elpd_loo), 
         ELPD_2 = sum(pw_2$elpd_loo),Delta_ELPD = sum(delta_i),SE_Delta = sqrt(length(delta_i) * var(delta_i)),
         Delta_over_SE = sum(delta_i) / sqrt(length(delta_i) * var(delta_i)),
         N_points_corrected = nrow(ex_1) + nrow(ex_2))
}

hybrid_full_reduced <- bind_rows(lapply(SECTOR_CODES, function(s)
  compare_hybrid(s, "Full_SORW", "Reduced_SORW")))
hybrid_sorw_forw <- bind_rows(lapply(SECTOR_CODES, function(s)
  compare_hybrid(s, "Reduced_FORW", "Reduced_SORW")))
hybrid_drift_sorw <- bind_rows(lapply(SECTOR_CODES, function(s)
  compare_hybrid(s, "Reduced_DriftFORW", "Reduced_SORW")))

classify_delta <- function(z) {
  case_when(z > 2 ~ "Favors Model 1", z < -2 ~ "Favors Model 2", TRUE ~ "Indistinguishable")
}

comparison_full_reduced <- step2_results$Full_vs_Reduced %>%
  transmute(Sector, PSIS_Delta = Delta_ELPD_Full_minus_Reduced,PSIS_SE = SE_Delta, PSIS_z = Delta_over_SE) %>%
  left_join(hybrid_full_reduced %>%transmute(Sector, Hybrid_Delta = Delta_ELPD, Hybrid_SE = SE_Delta,
                                             Hybrid_z = Delta_over_SE, N_points_corrected), by = "Sector") %>%
  mutate(PSIS_Conclusion = classify_delta(PSIS_z),Hybrid_Conclusion = classify_delta(Hybrid_z),
         Conclusion_Flipped = PSIS_Conclusion != Hybrid_Conclusion)

comparison_sorw_forw <- step2_results$SORW_vs_FORW %>%
  transmute(Sector, PSIS_Delta = Delta_ELPD_FORW_minus_SORW,PSIS_SE = SE_Delta, PSIS_z = Delta_over_SE) %>%
  left_join(hybrid_sorw_forw %>%transmute(Sector, Hybrid_Delta = Delta_ELPD, Hybrid_SE = SE_Delta,
                                          Hybrid_z = Delta_over_SE, N_points_corrected),by = "Sector") %>%
  mutate(PSIS_Conclusion = classify_delta(PSIS_z),Hybrid_Conclusion = classify_delta(Hybrid_z),
         Conclusion_Flipped = PSIS_Conclusion != Hybrid_Conclusion)

comparison_drift_sorw <- step2_results$Drift_vs_SORW %>%
  transmute(Sector, PSIS_Delta = Delta_ELPD_Drift_minus_SORW, PSIS_SE = SE_Delta, PSIS_z = Delta_over_SE) %>%
  left_join(hybrid_drift_sorw %>% transmute(Sector, Hybrid_Delta = Delta_ELPD, Hybrid_SE = SE_Delta,
                                            Hybrid_z = Delta_over_SE, N_points_corrected), by = "Sector") %>%
  mutate(PSIS_Conclusion = classify_delta(PSIS_z), Hybrid_Conclusion = classify_delta(Hybrid_z),
         Conclusion_Flipped = PSIS_Conclusion != Hybrid_Conclusion)

step8_exact_loo <- list(high_k_targets = high_k_targets, exact_results = exact_audit_results,
                        hybrid_full_reduced = hybrid_full_reduced, hybrid_sorw_forw = hybrid_sorw_forw,
                        hybrid_drift_sorw = hybrid_drift_sorw,comparison_full_reduced = comparison_full_reduced,
                        comparison_sorw_forw = comparison_sorw_forw,comparison_drift_sorw = comparison_drift_sorw)
saveRDS(step8_exact_loo, file.path(OUTPUT_DIR, "Step8_exact_LOO.rds"))
write_csv(comparison_drift_sorw, file.path(OUTPUT_DIR, "Step8_Comparison_Drift_vs_SORW.csv"))

cat("08_exact_LOO.R complete\n")
