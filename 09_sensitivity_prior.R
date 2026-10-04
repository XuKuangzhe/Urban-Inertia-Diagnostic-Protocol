# =============================================================================
# 09_sensitivity_prior.R
# Prior concentration / shrinkage sensitivity
# input (输入)：Step2_primary_results.rds
# output (输出)：Step9_prior_sensitivity.rds
# =============================================================================
message("=== 09_sensitivity_prior.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
prior_grid <- tribble(
  ~Prior_Setting,          ~prior_alpha, ~prior_beta, ~lasso_alpha, ~lasso_beta,
  "Tight (scale x0.5)",     1,            4,           30,           30,
  "Baseline (scale x1)",    1,            2,           30,           30,
  "Wide (scale x2)",        1,            1,           30,           30
)

extract_rpers <- function(fit) {
  post <- rstan::extract(fit, pars = c("s_mu", "s_Y"))
  R <- post$s_mu^2 / (post$s_mu^2 + post$s_Y^2)
  tibble(R_pers_Mean = mean(R), R_pers_Median = median(R),R_pers_Lower = quantile(R, 0.025), R_pers_Upper = quantile(R, 0.975))
}

analytic_p_exc_simple <- function(fit, Y_T, model = c("SORW", "FORW"),target_ratio = CFG$target_ratio, horizon = 12) {
  model <- match.arg(model)
  post <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y"))
  mu_trend <- post$mu_trend
  mu_T1 <- as.numeric(mu_trend[, ncol(mu_trend) - 1])
  mu_T  <- as.numeric(mu_trend[, ncol(mu_trend)])
  s_mu  <- as.numeric(post$s_mu)
  s_Y   <- as.numeric(post$s_Y)
  if (model == "FORW") {
    pred_mean <- mu_T
    pred_var  <- horizon * s_mu^2 + s_Y^2
  } else {
    pred_mean <- mu_T + horizon * (mu_T - mu_T1)
    pred_var  <- sum(seq_len(horizon)^2) * s_mu^2 + s_Y^2
  }
  threshold_log <- log(Y_T) + log(target_ratio)
  mean(1 - pnorm(threshold_log, mean = pred_mean, sd = sqrt(pred_var)))
}

model_comparison_rows  <- list()
persistence_rows       <- list()
forecast_rows          <- list()
diagnostic_rows        <- list()

total_tasks <- nrow(prior_grid) * length(SECTOR_CODES)
counter <- 0
cat(sprintf("\n>>> Starting Step 09: %d prior settings × %d sectors = %d total fits\n", 
            nrow(prior_grid), length(SECTOR_CODES), total_tasks))
pb <- utils::txtProgressBar(min = 0, max = total_tasks, style = 3)
for (p_idx in seq_len(nrow(prior_grid))) {
  prior_row <- prior_grid[p_idx, ]
  
  for (sector in SECTOR_CODES) {
    pretty <- SECTOR_LABELS[[sector]]
    Y <- as.numeric(D[[sector]])
    
    full_data    <- make_full_data(Y, prior = prior_row)
    reduced_data <- make_reduced_data(Y, prior = prior_row)
    seed_base <- 100000 + p_idx * 10000 + match(sector, SECTOR_CODES) * 100
    lab <- paste(prior_row$Prior_Setting, pretty)
    fit_full   <- fit_retry(MODELS$Full_SORW,         full_data,    seed_base + 1, paste(lab, "Full"))
    fit_rsorw  <- fit_retry(MODELS$Reduced_SORW,      reduced_data, seed_base + 2, paste(lab, "Reduced SORW"))
    fit_rforw  <- fit_retry(MODELS$Reduced_FORW,      reduced_data, seed_base + 3, paste(lab, "Reduced FORW"))
    fit_rdrift <- fit_retry(MODELS$Reduced_DriftFORW, reduced_data, seed_base + 4, paste(lab, "Reduced DriftFORW"))
    
    d_full   <- diagnose_stan_fit(fit_full,   "Full SORW",         pretty)
    d_rsorw  <- diagnose_stan_fit(fit_rsorw,  "Reduced SORW",      pretty)
    d_rforw  <- diagnose_stan_fit(fit_rforw,  "Reduced FORW",      pretty)
    d_rdrift <- diagnose_stan_fit(fit_rdrift, "Reduced DriftFORW", pretty)
    
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- bind_rows(
      d_full   %>% mutate(Prior_Setting = prior_row$Prior_Setting),
      d_rsorw  %>% mutate(Prior_Setting = prior_row$Prior_Setting),
      d_rforw  %>% mutate(Prior_Setting = prior_row$Prior_Setting),
      d_rdrift %>% mutate(Prior_Setting = prior_row$Prior_Setting))
    
    loo_full   <- loo(extract_log_lik(fit_full,   "log_lik", merge_chains = FALSE), cores = 1)
    loo_rsorw  <- loo(extract_log_lik(fit_rsorw,  "log_lik", merge_chains = FALSE), cores = 1)
    loo_rforw  <- loo(extract_log_lik(fit_rforw,  "log_lik", merge_chains = FALSE), cores = 1)
    loo_rdrift <- loo(extract_log_lik(fit_rdrift, "log_lik", merge_chains = FALSE), cores = 1)
    
    cmp_fr <- compare_pointwise(loo_full,   loo_rsorw, "Full SORW",         "Reduced SORW")
    cmp_rw <- compare_pointwise(loo_rforw,  loo_rsorw, "Reduced FORW",      "Reduced SORW")
    cmp_dr <- compare_pointwise(loo_rdrift, loo_rsorw, "Reduced DriftFORW", "Reduced SORW")
    
    model_comparison_rows[[length(model_comparison_rows) + 1]] <- bind_rows(
      tibble(Prior_Setting = prior_row$Prior_Setting, Sector = pretty,
             Comparison = "Full - Reduced", Delta_ELPD = cmp_fr$Delta_ELPD,
             SE_Delta = cmp_fr$SE_Delta, z = cmp_fr$Delta_over_SE),
      tibble(Prior_Setting = prior_row$Prior_Setting, Sector = pretty,
             Comparison = "FORW - SORW", Delta_ELPD = cmp_rw$Delta_ELPD,
             SE_Delta = cmp_rw$SE_Delta, z = cmp_rw$Delta_over_SE),
      tibble(Prior_Setting = prior_row$Prior_Setting, Sector = pretty,
             Comparison = "Drift - SORW", Delta_ELPD = cmp_dr$Delta_ELPD,
             SE_Delta = cmp_dr$SE_Delta, z = cmp_dr$Delta_over_SE))
    
    rpers <- extract_rpers(fit_rsorw)
    persistence_rows[[length(persistence_rows) + 1]] <- tibble(
      Prior_Setting = prior_row$Prior_Setting, Sector = pretty, !!!rpers)
    
    Y_2023 <- tail(Y, 1)
    p_sorw  <- analytic_p_exc_simple(fit_rsorw,  Y_2023, model = "SORW")
    p_forw  <- analytic_p_exc_simple(fit_rforw,  Y_2023, model = "FORW")
    p_drift <- analytic_p_exc(fit_rdrift, Y_2023, model = "DRIFT", horizons = 12)$P_exc
    
    forecast_rows[[length(forecast_rows) + 1]] <- tibble(
      Prior_Setting = prior_row$Prior_Setting, Sector = pretty,
      P_exc_SORW = p_sorw, P_exc_FORW = p_forw, P_exc_DRIFT = p_drift)
    
    rm(fit_full, fit_rsorw, fit_rforw, fit_rdrift, loo_full, loo_rsorw, loo_rforw, loo_rdrift); gc(verbose = FALSE)
    
    counter <- counter + 1
    utils::setTxtProgressBar(pb, counter)
  }
}
close(pb)

S9_model_comparison <- bind_rows(model_comparison_rows) %>%
  mutate(first  = sub(" - .*", "", Comparison), second = sub(".* - ", "", Comparison),
         Interpretation = case_when(z > 2 ~ paste0(first,  " favored"),z < -2 ~ paste0(second, " favored"),
                                    TRUE ~ "No clear predictive preference"))%>%select(-first, -second)

step9_prior <- list(prior_grid = prior_grid, model_comparison = S9_model_comparison,
                    persistence = bind_rows(persistence_rows), forecast = bind_rows(forecast_rows),
                    diagnostics = bind_rows(diagnostic_rows))
saveRDS(step9_prior, file.path(OUTPUT_DIR, "Step9_prior_sensitivity.rds"))

cat("09_sensitivity_prior.R complete\n")
# =============================================================================
# 09B_sensitivity_LASSO.R
# LASSO hyperprior sensitivity: only Full SORW depends on S_beta
# sigma prior fixed at CFG baseline; LASSO prior mean varied (not just concentration)
# =============================================================================
lasso_grid <- tribble(
  ~LASSO_Setting,                  ~lasso_alpha, ~lasso_beta,
  "Strong shrinkage (mean 0.25)",   30,           120,
  "Baseline (mean 1)",              30,            30,
  "Weak shrinkage (mean 2)",        30,            15
)

base_prior <- list(prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta,
                   lasso_alpha = CFG$lasso_alpha, lasso_beta = CFG$lasso_beta)
ctrl <- list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth)
cmp_rows <- list(); sbeta_rows <- list(); coef_rows <- list(); diag_rows <- list()
for (sector in SECTOR_CODES) {
  pretty <- SECTOR_LABELS[[sector]]
  Y <- as.numeric(D[[sector]])
  seed_base <- 200000 + match(sector, SECTOR_CODES) * 100
  # Reduced SORW reference (no LASSO term), fitted once per sector
  fit_rsorw <- fit_retry(MODELS$Reduced_SORW, make_reduced_data(Y, prior = base_prior),
                         seed_base + 1, paste(pretty, "Reduced SORW ref"))
  loo_rsorw <- loo(extract_log_lik(fit_rsorw, "log_lik", merge_chains = FALSE), cores = 1)
  diag_rows[[length(diag_rows) + 1]] <-
    diagnose_stan_fit(fit_rsorw, "Reduced SORW", pretty) %>% mutate(LASSO_Setting = "Reference")
  
  for (g in seq_len(nrow(lasso_grid))) {
    row <- lasso_grid[g, ]
    pr  <- modifyList(base_prior, list(lasso_alpha = row$lasso_alpha, lasso_beta = row$lasso_beta))
    fit_full <- fit_retry(MODELS$Full_SORW, make_full_data(Y, prior = pr),
                          seed_base + 10 + g, paste(pretty, row$LASSO_Setting))
    loo_full <- loo(extract_log_lik(fit_full, "log_lik", merge_chains = FALSE), cores = 1)
    k_full   <- loo_full$diagnostics$pareto_k
    diag_rows[[length(diag_rows) + 1]] <-
      diagnose_stan_fit(fit_full, "Full SORW", pretty) %>% mutate(LASSO_Setting = row$LASSO_Setting)
    # Full - Reduced (PSIS-LOO; exact-LOO replacement not repeated here)
    cmp <- compare_pointwise(loo_full, loo_rsorw, "Full SORW", "Reduced SORW")
    cmp_rows[[length(cmp_rows) + 1]] <- tibble(
      LASSO_Setting = row$LASSO_Setting, Sector = pretty,
      Delta_ELPD = cmp$Delta_ELPD, SE_Delta = cmp$SE_Delta, z = cmp$Delta_over_SE,
      N_pareto_k_gt_0_7 = sum(k_full > 0.7), Max_pareto_k = max(k_full))
    # S_beta posterior vs its prior
    sb <- as.numeric(rstan::extract(fit_full, pars = "S_beta")$S_beta)
    q  <- quantile(sb, c(0.025, 0.5, 0.975))
    prior_mean <- row$lasso_alpha / row$lasso_beta
    prior_sd   <- sqrt(row$lasso_alpha) / row$lasso_beta
    sbeta_rows[[length(sbeta_rows) + 1]] <- tibble(LASSO_Setting = row$LASSO_Setting, Sector = pretty,
                                                   Prior_Mean = prior_mean, Prior_SD = prior_sd,
                                                   Post_Mean = mean(sb), Post_Median = q[[2]], 
                                                   Post_Lower = q[[1]], Post_Upper = q[[3]],
                                                   Post_SD = sd(sb), SD_Contraction = 1 - sd(sb) / prior_sd)
    # Coefficients
    beta <- rstan::extract(fit_full, pars = "beta")$beta
    coef_rows[[length(coef_rows) + 1]] <- tibble(LASSO_Setting = row$LASSO_Setting, 
                                                 Sector = pretty, Driver = colnames(X_LOG),Mean = colMeans(beta),
                                                 Lower_95 = apply(beta, 2, quantile, 0.025), 
                                                 Upper_95 = apply(beta, 2, quantile, 0.975),
                                                 Prob_Direction = pmax(colMeans(beta > 0), colMeans(beta < 0))) %>%
      mutate(Excludes_Zero = Lower_95 > 0 | Upper_95 < 0)
    rm(fit_full, loo_full); gc()
  }
  rm(fit_rsorw, loo_rsorw); gc()
}

S9B_comparison <- bind_rows(cmp_rows) %>%
  mutate(Interpretation = case_when(z > 2 ~ "Full favored", z < -2 ~ "Reduced favored",
                                    TRUE ~ "No clear predictive preference"))
S9B_sbeta <- bind_rows(sbeta_rows)
S9B_coef  <- bind_rows(coef_rows)

# One-line summary per LASSO setting, for the SI table and Results text
S9B_coef_summary <- S9B_coef %>% group_by(LASSO_Setting) %>%
  summarise(N_coef = n(), N_CrI_excludes_zero = sum(Excludes_Zero),
            Max_Prob_Direction = max(Prob_Direction), Mean_abs_beta = mean(abs(Mean)), .groups = "drop")

step9B_lasso <- list(lasso_grid = lasso_grid, model_comparison = S9B_comparison,
                     S_beta = S9B_sbeta, coefficients = S9B_coef,
                     coef_summary = S9B_coef_summary, diagnostics = bind_rows(diag_rows))
saveRDS(step9B_lasso, file.path(OUTPUT_DIR, "Step9B_lasso_sensitivity.rds"))
write_csv(S9B_comparison,   file.path(OUTPUT_DIR, "Step9B_Full_vs_Reduced.csv"))
write_csv(S9B_sbeta,        file.path(OUTPUT_DIR, "Step9B_S_beta.csv"))
write_csv(S9B_coef_summary, file.path(OUTPUT_DIR, "Step9B_coef_summary.csv"))
cat("09B_sensitivity_LASSO.R complete\n")

cf <- step9B_lasso$coefficients %>% filter(LASSO_Setting == "Baseline (mean 1)")
sd_x <- purrr::map_dfr(SECTOR_CODES, function(s) {
  X <- step2_data_v3[[s]]$Full_SORW$X          # 若 X 的元素名不同请改
  tibble(Sector = SECTOR_LABELS[[s]], Driver = c("P", "A", "Tstr", "Ttech"), SD_X = apply(X, 2, sd))
})
sy <- sigma_pp %>% select(Sector, s_Y_post_med)
mde <- cf %>% left_join(sd_x, by = c("Sector", "Driver")) %>% left_join(sy, by = "Sector") %>%
  mutate(Beta_bound = pmax(abs(Lower_95), abs(Upper_95)),
         Max_effect_1SD = Beta_bound * SD_X,
         Ratio_to_sigmaY = Max_effect_1SD / s_Y_post_med)
print(as.data.frame(mde))