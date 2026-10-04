# =============================================================================
# 03_identification_diagnostics.R
# 03A Posterior R_pers (R_pers 后验)
# 03B R_pers prior -> posterior
# 03C beta posterior(beta 后验)
# 03D S_beta posterior(S_beta 后验)
# 03E Full vs Reduced predictive evidence (Full vs Reduced 预测证据)
# 输入：Step2_primary_results.rds
# 输出：Step3_identification.rds
# =============================================================================
message("=== 03_identification_diagnostics.R ===")
#step2_results <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
# -----------------------------------------------------------------------------
# 03A Posterior R_pers (R_pers 后验)
# -----------------------------------------------------------------------------
persistence_rows <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  fit_sorw <- step2_results$fits[[sector]]$Reduced_SORW
  post <- rstan::extract(fit_sorw, pars = c("s_mu", "s_Y"))
  sigma_mu <- as.numeric(post$s_mu)
  sigma_Y  <- as.numeric(post$s_Y)
  stopifnot(all(sigma_mu > 0), all(sigma_Y > 0))
  q<- sigma_mu^2 / sigma_Y^2
  R_pers <- sigma_mu^2 / (sigma_mu^2 + sigma_Y^2)
  mu_s <- summarize_posterior(sigma_mu)
  Y_s  <- summarize_posterior(sigma_Y)
  q_s  <- summarize_posterior(q)
  R_s  <- summarize_posterior(R_pers)
  
  rw_row <- step2_results$SORW_vs_FORW %>% filter(Sector == SECTOR_LABELS[[sector]])
  stopifnot(nrow(rw_row) == 1)
  
  tibble(Sector = SECTOR_LABELS[[sector]],Sigma_mu_Mean = mu_s$Mean, Sigma_mu_Median = mu_s$Median,
         Sigma_mu_Lower = mu_s$Lower, Sigma_mu_Upper = mu_s$Upper,
         Sigma_Y_Mean = Y_s$Mean, Sigma_Y_Median = Y_s$Median,
         Sigma_Y_Lower = Y_s$Lower, Sigma_Y_Upper = Y_s$Upper,
         q_Mean = q_s$Mean, q_Median = q_s$Median,
         q_Lower = q_s$Lower, q_Upper = q_s$Upper,
         R_pers_Mean = R_s$Mean, R_pers_Median = R_s$Median,
         R_pers_Lower = R_s$Lower, R_pers_Upper = R_s$Upper,
         P_q_gt_1 = mean(q > 1),
         P_R_pers_gt_0_5 = mean(R_pers > 0.5),
         Delta_ELPD_FORW_minus_SORW = rw_row$Delta_ELPD_FORW_minus_SORW,
         SE_Delta = rw_row$SE_Delta,
         Delta_over_SE = rw_row$Delta_over_SE,
         Evidence = rw_row$Evidence
  )
})

persistence_table<-persistence_rows%>%mutate(Sector=factor(Sector,levels=SECTOR_ORDER_LABELS))%>%arrange(Sector)

# -----------------------------------------------------------------------------
# 03B Prior-posterior R_pers
# -----------------------------------------------------------------------------
set.seed(20260924)
N_PRIOR <- 200000
prior_sigma_mu <- rlnorm(N_PRIOR, meanlog = -2.0, sdlog = 0.7)
prior_sigma_Y  <- rlnorm(N_PRIOR, meanlog = -2.3, sdlog = 0.6)
prior_Rpers    <- prior_sigma_mu^2 / (prior_sigma_mu^2 + prior_sigma_Y^2)
prior_q        <- quantile(prior_Rpers, probs = c(0.025, 0.5, 0.975))

posterior_Rpers <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  fit  <- step2_results$fits[[sector]]$Reduced_SORW
  post <- rstan::extract(fit, pars = c("s_mu", "s_Y"))
  tibble(Sector = SECTOR_LABELS[[sector]],R_pers = post$s_mu^2 / (post$s_mu^2 + post$s_Y^2))
})

Rpers_prior_posterior <- posterior_Rpers %>%group_by(Sector) %>%
  summarise(Posterior_Median = median(R_pers), Posterior_Lower = quantile(R_pers, 0.025),
            Posterior_Upper = quantile(R_pers, 0.975), Posterior_SD = sd(R_pers), Prior_Median = prior_q[2], 
            Prior_Lower = prior_q[1], Prior_Upper = prior_q[3], Prior_SD = sd(prior_Rpers),
            CI_Width_Contraction = 1 - ((Posterior_Upper - Posterior_Lower) / (prior_q[3] - prior_q[1])),
            Median_Shift = Posterior_Median - prior_q[2], .groups = "drop")

# -----------------------------------------------------------------------------
# 03C beta posterior (beta 后验)
# -----------------------------------------------------------------------------
driver_names <- c("Population", "Affluence", "Structure", "Technology")

coef_rows <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  fit_full <- step2_results$fits[[sector]]$Full_SORW
  post <- rstan::extract(fit_full, pars = c("beta", "beta_covid", "S_beta"))
  beta_draws <- post$beta
  stopifnot(ncol(beta_draws) == 4)
  
  rows <- purrr::map_dfr(1:4, function(j) {
    x <- beta_draws[, j]
    q <- quantile(x, probs = c(0.025, 0.5, 0.975))
    p_pos <- mean(x > 0)
    tibble(Sector = SECTOR_LABELS[[sector]],Driver = driver_names[j],Mean = mean(x),Median = q[2],
           Lower_95 = q[1],Upper_95 = q[3],Prob_Positive = p_pos,Prob_Negative = 1 - p_pos,
           Prob_Direction = max(p_pos, 1 - p_pos), Zero_in_95_CrI = q[1] <= 0 & q[3] >= 0)
  })
  
  covid_x <- post$beta_covid
  cq <- quantile(covid_x, probs = c(0.025, 0.5, 0.975))
  cp <- mean(covid_x > 0)
  attr(rows, "covid") <- tibble(Sector = SECTOR_LABELS[[sector]],
                                Mean = mean(covid_x), Lower_95 = cq[1], Upper_95 = cq[3],
                                Prob_Direction = max(cp, 1 - cp), Zero_in_95_CrI = cq[1] <= 0 & cq[3] >= 0)
  rows
})

coef_covid <- coef_rows %>% group_split(Sector) %>% purrr::map_dfr(~ attr(.x, "covid"))
coef_table <- coef_rows %>% mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS),
                                   Driver = factor(Driver, levels = driver_names)) %>% arrange(Sector, Driver)

# -----------------------------------------------------------------------------
# 03D S_beta posterior (S_beta 后验)
# -----------------------------------------------------------------------------
S_beta_table <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  fit_full <- step2_results$fits[[sector]]$Full_SORW
  post <- rstan::extract(fit_full, pars = "S_beta")
  x <- as.numeric(post$S_beta)
  q <- quantile(x, probs = c(0.025, 0.5, 0.975))
  tibble(Sector = SECTOR_LABELS[[sector]], S_beta_Mean = mean(x), S_beta_Median = q[2],
         S_beta_Lower = q[1], S_beta_Upper = q[3], S_beta_SD = sd(x))
})

# -----------------------------------------------------------------------------
# 03E Full vs Reduced （预测证据）
# -----------------------------------------------------------------------------
loo_covariate_table<-step2_results$Full_vs_Reduced %>%mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS))

# =============================================================================
# Save (保存)
# =============================================================================
step3_results <- list(
  metadata = list(persistence_definition = "R_pers = sigma_mu^2 / (sigma_mu^2 + sigma_Y^2)",
                  q_definition = "q = sigma_mu^2 / sigma_Y^2",
                  interpretation = paste("Temporal noise-ratio diagnostic under fitted Reduced SORW.",
                                         "Not a causal or physical lock-in index.")),
  persistence_table= persistence_table,Rpers_prior_posterior= Rpers_prior_posterior,
  prior_Rpers= prior_Rpers,posterior_Rpers= posterior_Rpers,coefficient_table= coef_table,
  coefficient_covid= coef_covid,S_beta_table= S_beta_table,loo_covariate_table= loo_covariate_table
)

saveRDS(step3_results, file.path(OUTPUT_DIR, "Step3_identification.rds"))
write_csv(persistence_table,     file.path(OUTPUT_DIR, "Step3_Temporal_Persistence.csv"))
write_csv(coef_table,            file.path(OUTPUT_DIR, "Step4_Socioeconomic_Coefficients.csv"))
write_csv(S_beta_table,          file.path(OUTPUT_DIR, "Step4_S_beta_Posterior.csv"))

cat("03_identification_diagnostics.R complete\n")
