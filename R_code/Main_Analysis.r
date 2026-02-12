# ==============================================================================
# Script Name: Main_Analysis.R
# Description: Generates main figures (Fig 1-5) for TFSC submission.
#              Includes Model Comparison, Sensitivity Analysis, and Projections.
# Author: XU Kuangzhe
# ==============================================================================

# ==============================================================================
# Part 0: Environment Setup & Data Loading
# ==============================================================================

library(rstan)
library(loo)        # For LOO-CV model comparison
library(tidyverse)  # For data manipulation and plotting
library(tidybayes)  # For extracting Bayesian credible intervals
source("~/HDI.R")

# Parallel computing settings
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# --- 1. Load Compiled Stan Models ---
# Ensure paths point to your local .stan files
NullSSM <- stan_model("~/PATH/SSM_Null_Dynamic.stan")
LSWaicM <- stan_model("~/PATH/SSM_Full_Dynamic.stan")

# --- 2. Data Preparation ---
# d: 能源消耗矩阵 (Rows: Years, Cols: Sectors)
# dipat: STIRPAT 驱动因子 (Population, GDP, Structure, Tech)
# PowN: 能源名称向量
dipat<-read_csv("~/PATH/tableA4.csv")
d <- read.csv("~/PATH/tableA5.csv")
# Create COVID Dummy Vector (Assuming T=14, 2010-2023)
# 2020, 2021, 2022 are typically considered impacted years in China
# Adjust this based on your specific sector data behavior (e.g., just 2020 for Oil?)


# Define Intervention Variable (COVID-19 Shock)
# Assumption: Pandemic impact covers 2020, 2021, and 2022.
if(!exists("covid_vec")) {
  covid_vec <- rep(0, 14) # Assuming N=14 (2010-2023)
  names(covid_vec) <- 2010:2023
  covid_vec[names(covid_vec) %in% c("2020", "2021", "2022")] <- 1
}

# Mapping for Professional Sector Labels in Plots
energy_labels <- list(
  "coal"     = "Coal",
  "oil"      = "Oil Products",
  "NGas"     = "Natural Gas",
  "LiNGas"   = "LPG",
  "thermal"  = "Thermal (Heat)",
  "electric" = "Electricity",
  "other"    = "Other"
)

# ==============================================================================
# Part 1: Table 1 - Monte Carlo Parameter Recovery 
# ==============================================================================

# 1. Estimate Ground Truth parameters from empirical data
# ------------------------------------------------------------------------------
cat("Estimating Ground Truth parameters from empirical data...\n")

sumTrue <- data.frame()
# Assuming columns 2 to 8 represent energy sectors
energy_cols <- colnames(d)[2:8] 

for(energy_name in energy_cols){
  # Fit model to empirical data to get posterior means
  datalist_null <- list(T=14, Y=d[[energy_name]], covid_dummy=covid_vec, prior_alpha=20, prior_beta=20) 
  fit_temp <- sampling(NullSSM, data=datalist_null, iter=4000, seed=1234, refresh=0)
  
  post <- rstan::extract(fit_temp)
  dture <- tibble(
    Energy = energy_name,
    true_s_Y = mean(post$s_Y),
    true_s_mu = mean(post$s_mu),
    true_beta = mean(post$beta_covid)
  )
  sumTrue <- bind_rows(sumTrue, dture)
}

param_df <- sumTrue
print("Ground Truth Parameters defined:")
print(param_df)

# 2. Run Simulation Experiment
# ------------------------------------------------------------------------------
T_sim <- 14     # Time series length N=14
N_sims <- 100   # Number of simulations (Recommended: 100)
final_results <- data.frame()

for(i in 1:nrow(param_df)) {
  
  curr_energy <- param_df$Energy[i]
  t_s_mu <- param_df$true_s_mu[i]
  t_s_Y  <- param_df$true_s_Y[i]
  t_beta <- param_df$true_beta[i]
  
  cat(paste0("\nSimulating: ", curr_energy, " (s_mu=", round(t_s_mu,3), ")\n"))
  
  temp_coverage <- data.frame()
  
  for(sim in 1:N_sims) {
    # --- A. Generate Synthetic Data ---
    set.seed(sim * 100 + i)
    
    mu_trend <- numeric(T_sim)
    mu_trend[1] <- 5.0 
    mu_trend[2] <- rnorm(1, mean = mu_trend[1], sd = t_s_mu)
    
    # 2nd-order Random Walk for trend
    for(t in 3:T_sim) {
      mu_trend[t] <- 2*mu_trend[t-1] - mu_trend[t-2] + rnorm(1, 0, t_s_mu)
    }
    
    # Add Pandemic Shock
    mu_total <- mu_trend + t_beta * covid_vec
    
    # Generate Observations (log scale)
    y_log <- rnorm(T_sim, mean = mu_total, sd = t_s_Y)
    
    # --- B. Refit Stan Model ---
    # Pass exp(y_log) as the Stan model takes raw values and logs internally
    datalist <- list(T=T_sim, Y=exp(y_log), covid_dummy=covid_vec, prior_alpha=20, prior_beta=20)
    
    # Fast sampling for validation (iter=2000)
    fit_sim <- sampling(NullSSM, data=datalist, iter=2000, refresh=0)
    post_sim <- rstan::extract(fit_sim)
    
    # --- C. Check Coverage (95% CI) ---
    # Check s_mu (Trend Noise)
    ci_mu <- quantile(post_sim$s_mu, c(0.025, 0.975))
    cover_mu <- (t_s_mu >= ci_mu[1] & t_s_mu <= ci_mu[2])
    
    # Check s_Y (Obs Noise)
    ci_Y <- quantile(post_sim$s_Y, c(0.025, 0.975))
    cover_Y <- (t_s_Y >= ci_Y[1] & t_s_Y <= ci_Y[2])
    
    # Check beta_covid (Shock Magnitude)
    ci_beta <- quantile(post_sim$beta_covid, c(0.025, 0.975))
    cover_beta <- (t_beta >= ci_beta[1] & t_beta <= ci_beta[2])
    
    temp_coverage <- bind_rows(temp_coverage, 
                               data.frame(Sim = sim, 
                                          Cover_Mu = cover_mu, 
                                          Cover_Y = cover_Y,
                                          Cover_Beta = cover_beta))
    
    if(sim %% 20 == 0) cat(".")
  }
  
  # Aggregate results for current sector
  stats <- temp_coverage %>% 
    summarise(
      Energy = curr_energy,
      Coverage_Trend_Noise = mean(Cover_Mu),
      Coverage_Obs_Noise   = mean(Cover_Y),
      Coverage_Shock_Beta  = mean(Cover_Beta)
    )
  final_results <- bind_rows(final_results, stats)
}

print("=========== Table 1: Parameter Recovery Results ===========")
print(final_results)
# write_csv(final_results, "Table1_ParameterRecovery.csv")


# ==============================================================================
# Part 2: Figure 1 - Sensitivity to Time Series Length 
# ==============================================================================
# 1. Simulation Setup
# ------------------------------------------------------------------------------
N_levels <- c(10, 12, 14, 16, 20, 25) # Tested sample sizes
N_sims_per_level <- 100               # Simulations per level

# Set Ground Truth (Using Oil Sector parameters as a high-volatility stress test)
true_pars <- list(
  s_mu = 0.734,   # Trend Noise
  s_Y  = 0.642,   # Obs Noise
  beta = -0.155   # Shock Magnitude
)

sensitivity_results <- data.frame()

# 2. Loop through different lengths N
# ------------------------------------------------------------------------------
for(n_curr in N_levels) {
  
  cat(paste0("Testing Length N = ", n_curr, "... "))
  
  # Dynamic Covid Vector (Ensure shock is always at the end, e.g., last 3 years)
  covid_vec_sim <- rep(0, n_curr)
  impact_idx <- c(n_curr-3, n_curr-2, n_curr-1)
  impact_idx <- impact_idx[impact_idx > 0]
  covid_vec_sim[impact_idx] <- 1
  
  temp_res <- data.frame()
  
  for(sim in 1:N_sims_per_level) {
    set.seed(sim * 1000 + n_curr)
    
    # Generate Synthetic Data
    mu_trend <- numeric(n_curr)
    mu_trend[1] <- 5.0
    mu_trend[2] <- rnorm(1, mu_trend[1], true_pars$s_mu)
    for(t in 3:n_curr) mu_trend[t] <- 2*mu_trend[t-1] - mu_trend[t-2] + rnorm(1, 0, true_pars$s_mu)
    
    mu_total <- mu_trend + true_pars$beta * covid_vec_sim
    y_log <- rnorm(n_curr, mu_total, true_pars$s_Y)
    
    # Fit Model
    d_sim <- list(T=n_curr, Y=exp(y_log), covid_dummy=covid_vec_sim, prior_alpha=20, prior_beta=20)
    fit_sim <- sampling(NullSSM, data=d_sim, iter=2000, chains=2, refresh=0)
    post <- rstan::extract(fit_sim)
    
    # Check Coverage
    ci_beta <- quantile(post$beta_covid, c(0.025, 0.975))
    cover_beta <- (true_pars$beta >= ci_beta[1] & true_pars$beta <= ci_beta[2])
    
    ci_mu <- quantile(post$s_mu, c(0.025, 0.975))
    cover_mu <- (true_pars$s_mu >= ci_mu[1] & true_pars$s_mu <= ci_mu[2])
    
    temp_res <- bind_rows(temp_res, data.frame(Cover_Beta=cover_beta, Cover_Mu=cover_mu))
  }
  
  # Aggregate Stats
  stats <- temp_res %>% summarise(
    N = n_curr,
    Coverage_Beta = mean(Cover_Beta),
    Coverage_Trend = mean(Cover_Mu)
  )
  sensitivity_results <- bind_rows(sensitivity_results, stats)
  cat("Done.\n")
}

# 3. Plot Sensitivity Results
# ------------------------------------------------------------------------------
p_sens <- ggplot(sensitivity_results, aes(x = N)) +
  # Shock Parameter Coverage (Solid Line)
  geom_line(aes(y = Coverage_Beta, color = "Shock Parameter"), linewidth = 1.2) +
  geom_point(aes(y = Coverage_Beta, color = "Shock Parameter"), size = 3) +
  
  # Trend Noise Coverage (Dashed Line)
  geom_line(aes(y = Coverage_Trend, color = "Trend Noise"), linewidth = 1.2, linetype="dashed") +
  geom_point(aes(y = Coverage_Trend, color = "Trend Noise"), size = 3) +
  
  # 95% Reference Line
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "gray30") +
  
  scale_color_manual(values = c("Shock Parameter" = "#D55E00", "Trend Noise" = "#0072B2")) +
  scale_x_continuous(breaks = N_levels) +
  ylim(0.8, 1.05) +
  
  labs(title = "Sensitivity of Parameter Recovery to Time Series Length",
       x = "Time Series Length (N)",
       y = "Coverage Probability (95% HDI)",
       color = "Parameter") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(p_sens)
# ggsave("Figure 1_LengthSensitivity.png", plot=p_sens, width=8, height=5)

# ==============================================================================
# Part 3: Figure 2 - Model Comparison (Full vs. Null)
# ==============================================================================

loo_comparison_list <- list()  
loo_summary_df <- data.frame() 

# Loop through energy sectors (Assuming cols 2 to 8 contain sector data)
for(o in 2:8){
  energy_name <- colnames(d)[o] 
  print(paste("Processing LOO comparison for:", energy_name))
  
  # --- Data Pre-processing ---
  # Log-transformation of STIRPAT predictors
  X_log <- dipat[, 2:5] %>% mutate(covid_dummy = covid_vec)
  X_log$P <- log(X_log$P)
  X_log$A <- log(X_log$A)
  X_log$Tstr <- log(X_log$Tstr)
  X_log$Ttech <- log(X_log$Ttech)
  
  # --- Prepare Stan Data Lists ---
  # Full Model: Includes socioeconomic drivers X
  datalist_full <- list(T=14, K=5, X=X_log, Y=d[,o], prior_alpha=20, prior_beta=20)
  # Null Model: Inertia only + Intervention
  datalist_null <- list(T=14, Y=d[,o], covid_dummy=covid_vec, prior_alpha=20, prior_beta=20) 
  
  # --- Sampling ---
  fit1 <- sampling(LSWaicM, data=datalist_full, iter=4000, seed=1234, refresh=0)
  fit2 <- sampling(NullSSM, data=datalist_null, iter=4000, seed=1234, refresh=0)
  
  # --- LOO Calculation ---
  log_lik_full <- extract_log_lik(fit1, merge_chains = FALSE)
  log_lik_null <- extract_log_lik(fit2, merge_chains = FALSE)
  
  loo_full <- loo(log_lik_full)
  loo_null <- loo(log_lik_null)
  
  # --- Comparison ---
  comp <- loo_compare(list(Full_Model = loo_full, Null_Model = loo_null))
  
  # --- Result Aggregation ---
  comp_df <- as.data.frame(comp) %>%
    rownames_to_column(var = "Model_Name") %>%
    mutate(Energy = energy_name)
  
  loo_summary_df <- rbind(loo_summary_df, comp_df)
}

# --- Plotting Figure 1 ---
plot_data_f1 <- loo_summary_df %>%
  filter(Model_Name == "Full_Model") %>% # Extract difference relative to Null Model
  mutate(
    Energy_Label = recode(Energy, !!!energy_labels),
    Energy_Label = reorder(Energy_Label, elpd_diff),
    lower = elpd_diff - 2 * se_diff,
    upper = elpd_diff + 2 * se_diff,
    is_significant = upper < 0 # Check if significantly worse than Null Model
  )

p1 <- ggplot(plot_data_f1, aes(x = elpd_diff, y = Energy_Label)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_col(aes(fill = is_significant), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#99C9E7"), guide = "none") +
  labs(
    title = "Model Comparison (LOO-CV)",
    subtitle = "Predictive difference: Full Model (Growth) vs. Null Model (Inertia)",
    x = expression(paste(Delta, "ELPD (Full Model - Null Model)")),
    y = ""
  ) +
  annotate("text", x = -0.5, y = 0.5, label = "Favors Null Model\n(Inertia)", 
           hjust = 1, size = 3.5, fontface = "italic", color = "gray30") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major.y = element_blank())

print(p1)
# ggsave("Figure 2_ModelComparison.png", plot=p1, width=9, height=5)


# ==============================================================================
# Part 4: Figure 3 - Sensitivity Analysis (Prior Specification)
# ==============================================================================

# Define Prior Scenarios
Prior_Scenarios <- list(
  "Weak"     = list(a = 0.01, b = 0.01, label = "Weak Priors (Gamma 0.01, 0.01)"),
  "Moderate" = list(a = 2.0,  b = 2.0,  label = "Moderate Priors (Gamma 2, 2)"),
  "Strong"   = list(a = 20.0, b = 20.0, label = "Strong Priors (Gamma 20, 20)")
)

results_list <- list()

for (energy in colnames(d)[2:8]) {
  cat(paste0("\n=== Processing Sensitivity: ", energy, " ===\n"))
  
  # Data Setup
  X_curr <- dipat[, 2:5] %>% mutate(covid_dummy = covid_vec)
  X_curr[,1:4] <- log(X_curr[,1:4])
  Y_curr <- d[[energy]]
  
  for (scenario_name in names(Prior_Scenarios)) {
    scen <- Prior_Scenarios[[scenario_name]]
    
    # Construct Data Lists with varying priors
    d_full <- list(T=14, K=5, X=X_curr, Y=Y_curr, prior_alpha=scen$a, prior_beta=scen$b)
    d_null <- list(T=14, Y=Y_curr, covid_dummy=covid_vec, prior_alpha=scen$a, prior_beta=scen$b)
    
    # Fit Models
    fit_full <- sampling(LSWaicM, data = d_full, iter = 4000, seed = 1234, refresh = 0)
    fit_null <- sampling(NullSSM, data = d_null, iter = 4000, seed = 1234, refresh = 0)
    
    # Calculate LOO difference point-wise to get SE
    loo_full <- loo(fit_full)
    loo_null <- loo(fit_null)
    elpd_diff_pt <- loo_full$pointwise[,"elpd_loo"] - loo_null$pointwise[,"elpd_loo"]
    
    results_list[[length(results_list) + 1]] <- tibble(
      Sector = energy,
      Prior_Strength = scen$label,
      Delta_ELPD = sum(elpd_diff_pt),
      SE = sqrt(length(elpd_diff_pt) * var(elpd_diff_pt))
    )
  }
}

# --- Plotting Figure 2 ---
plot_data_f2 <- bind_rows(results_list) %>%
  mutate(
    Sector = recode(Sector, !!!energy_labels),
    Prior_Strength = factor(Prior_Strength, levels = c(
      "Weak Priors (Gamma 0.01, 0.01)", 
      "Moderate Priors (Gamma 2, 2)", 
      "Strong Priors (Gamma 20, 20)"
    ))
  )

p2 <- ggplot(plot_data_f2, aes(x = Delta_ELPD, y = Sector, color = Prior_Strength)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.6) +
  geom_errorbarh(aes(xmin = Delta_ELPD - 2*SE, xmax = Delta_ELPD + 2*SE), 
                 height = 0.3, position = position_dodge(width = 0.7)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  labs(
    title = "Sensitivity Analysis",
    x = expression(paste(Delta, "ELPD (Full Model - Null Model)")), y = "",
    color = "Regularization Strength"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.border = element_rect(color="black", fill=NA))

print(p2)
# ggsave("Figure4_Sensitivity.png", plot=p2, width=11, height=7)


# ==============================================================================
# Part 5: Figures 4, 5, 6 - Structural Analysis & Projections
# ==============================================================================

output_dir <- "~/PATH/Optimized_Plots/"
if(!dir.exists(output_dir)) dir.create(output_dir)

for(o in 2:8){
  
  raw_name <- colnames(d)[o]
  energy_name <- ifelse(!is.null(energy_labels[[raw_name]]), energy_labels[[raw_name]], raw_name)
  print(paste(">>> Analyzing Structural Dynamics:", energy_name))
  
  # --- 1. Fit Best Model (Null Model) ---
  datalist_null <- list(T=14, Y=d[,o], covid_dummy=covid_vec, prior_alpha=20, prior_beta=20)
  fit_obj <- sampling(NullSSM, data=datalist_null, iter=4000, seed=1234, refresh=0)
  
  # Extract Parameters
  ms <- rstan::extract(fit_obj)
  y_obs <- log(d[,o])
  T_obs <- length(y_obs)
  
  # --- 2. Figure 4: Decoupling Index (DI) Evolution ---
  # Calculate Rolling Signal-to-Noise Ratio
  # Remove Pandemic Effect to see underlying inertia
  impact_matrix <- sweep(matrix(rep(covid_vec, each=4000), nrow=4000), 1, ms$beta_covid, "*")
  fitted_trend <- colMeans(ms$mu_trend)
  fitted_noise <- y_obs - (fitted_trend + colMeans(impact_matrix)) 
  
  di_seq <- c()
  window_size <- 3
  for(i in 1:T_obs){
    idx <- max(1, i - window_size + 1):i
    ss_trend <- sum((fitted_trend[idx] - mean(fitted_trend[idx]))^2)
    ss_noise <- sum(fitted_noise[idx]^2)
    di_seq <- c(di_seq, ss_trend / (ss_trend + ss_noise + 1e-6))
  }
  
  p3 <- ggplot(data.frame(Year=2010:2023, DI=di_seq), aes(x=Year, y=DI)) +
    geom_smooth(method = "loess", se = FALSE, color = "#08519c", alpha = 0.2, linewidth = 2) +
    geom_line(color = "#08519c", linewidth = 1) +
    geom_point(size = 3, color = "#08519c", fill = "white", shape = 21, stroke = 1.5) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.1)) +
    labs(title = paste0("Inertia Evolution: ", energy_name), y = "Decoupling Index (DI)", x = "") +
    theme_classic(base_size = 14)
  
  ggsave(paste0(output_dir, "Fig5_DI_", raw_name, ".png"), plot = p3, width = 6, height = 4)
  
  # --- 3. Generate Future Projections ---
  pred_years <- 12 # 2024-2035
  n_iter <- length(ms$s_mu)
  y_future <- matrix(NA, nrow = n_iter, ncol = pred_years)
  mu_curr <- ms$mu_trend[, T_obs]
  mu_prev <- ms$mu_trend[, T_obs - 1]
  
  set.seed(1234)
  for (t in 1:pred_years) {
    # 2nd-order Random Walk extrapolation
    mu_next <- rnorm(n_iter, mean = 2 * mu_curr - mu_prev, sd = ms$s_mu)
    y_future[, t] <- rnorm(n_iter, mean = mu_next, sd = ms$s_Y)
    mu_prev <- mu_curr
    mu_curr <- mu_next
  }
  
  # --- 5. Figure 4: Probabilistic Fan Chart ---
  # Combine Historical and Forecast Data
  df_fan <- bind_rows(
    as.data.frame(ms$mu_trend) %>% mutate(Type="Hist", Year_Start=2010),
    as.data.frame(y_future) %>% mutate(Type="Pred", Year_Start=2024)
  ) %>%
    pivot_longer(cols = starts_with("V"), names_to = "idx", values_to = "val") %>%
    mutate(Year = Year_Start + as.numeric(gsub("V","",idx)) - 1) %>%
    group_by(Year, Type) %>%
    median_qi(val, .width = c(0.5, 0.8, 0.95))
  
  p4 <- ggplot(df_fan, aes(x = Year, y = val)) +
    geom_ribbon(data=filter(df_fan, Type=="Pred"), aes(ymin=val.lower, ymax=val.upper, fill=as.factor(.width)), alpha=0.8) +
    geom_line(data=filter(df_fan, Type=="Hist"), aes(y=val), color="black") +
    geom_line(data=filter(df_fan, Type=="Pred"), aes(y=val), color="#08519c") +
    geom_point(data=data.frame(Year=2010:2023, val=y_obs), shape=21, fill="white") +
    scale_fill_brewer(palette = "Blues") +
    labs(title = paste0("Projection: ", energy_name), y = "Log Consumption", x = "") +
    theme_classic(base_size = 14) + theme(legend.position="none")
  
  ggsave(paste0(output_dir, "Fig6_Fan_", raw_name, ".png"), plot = p4, width = 6, height = 4)
  
  # --- 5. Figure 6: Risk Density Plot ---
  # Define Policy Target (e.g., 20% reduction by 2035)
  target_val <- y_obs[length(y_obs)] + log(0.8) 
  y_2035 <- y_future[, 12]
  prob_fail <- mean(y_2035 > target_val)
  
  dens <- density(y_2035, adjust = 1.5)
  df_dens <- data.frame(x = dens$x, y = dens$y)
  
  p5 <- ggplot(df_dens, aes(x, y)) +
    geom_area(fill = "grey90") +
    geom_area(data = filter(df_dens, x > target_val), fill = "#de2d26", alpha = 0.8) +
    geom_vline(xintercept = target_val, linetype = "dashed") +
    annotate("text", x = target_val, y = max(dens$y), label = "Target", vjust = -0.5) +
    labs(title = paste0("Risk Analysis: ", energy_name), 
         subtitle = paste0("Failure Probability: ", round(prob_fail*100, 1), "%"),
         x = "2035 Log Consumption", y = "Density") +
    theme_classic(base_size = 14)
  
  ggsave(paste0(output_dir, "Fig7_Risk_", raw_name, ".png"), plot = p5, width = 6, height = 4)
}

# ==============================================================================
# Part 6: Figures 7 - Structural Intervention "Stress Test" (Multi-Scenario)
# ==============================================================================
# Objective: Test the impact of varying intensities of structural intervention 
# (0%, 10%, 20%, 30%, 40%) on the Probability of Target Exceedance in 2035.
# ------------------------------------------------------------------------------

library(ggridges) # Required for Ridge Plots

target_sector_idx <- which(colnames(d) == "oil") 
energy_name_risk <- "Oil Products"

# 1. Data and Model Preparation
# ------------------------------------------------------------------------------
# Reuse existing fit object if available, otherwise re-sample.
if(!exists("ms_risk")) {
  datalist_risk <- list(T=14, Y=d[,target_sector_idx], covid_dummy=covid_vec, prior_alpha=20, prior_beta=20)
  fit_risk <- sampling(NullSSM, data=datalist_risk, iter=4000, seed=1234, refresh=0)
  ms_risk <- rstan::extract(fit_risk)
}

# 2. Define Scenario Parameters
# ------------------------------------------------------------------------------
pred_years <- 12 # Projection horizon: 2024-2035
year_start <- 2024
intervention_year <- 2025

# Define Intervention Intensities (Shock Levels)
# From 0% (BAU) to 40% (Deep Exnovation)
shock_levels <- c(0, 0.10, 0.20, 0.30, 0.40)
shock_labels <- paste0("-", shock_levels * 100, "% (Regime Shift)")
shock_labels[1] <- "BAU (High Inertia)"

# Define Regime Stability Factors (Volatility Damping)
# Logic: Mild interventions (-10%) do not alter the regime (low damping).
# Deep interventions (-40%) imply a shift to a regulated utility model (high damping, e.g., 90%).
# This models a non-linear stability gain.
vol_damp_factors <- c(0, 0.1, 0.5, 0.8, 0.9) 

# Define Policy Cap (Target): 20% reduction from 2023 levels
target_cap <- log(d[14, target_sector_idx]) + log(0.80) 

# 3. Run Simulation Loop
# ------------------------------------------------------------------------------
summary_df <- data.frame()
results_df <- data.frame()

set.seed(1234) # Ensure reproducibility

for (i in 1:length(shock_levels)) {
  lvl <- shock_levels[i]
  damp <- vol_damp_factors[i]
  lbl <- shock_labels[i]
  magnitude <- log(1 - lvl)
  
  # --- Core Mechanism: Regime Switch ---
  # Under deep intervention, the system shifts from a stochastic random walk 
  # to a planned trajectory, compressing the posterior variance.
  reduced_s_mu <- ms_risk$s_mu * (1 - damp)
  reduced_s_Y  <- ms_risk$s_Y  * (1 - damp)
  
  cat(paste0("Simulating: ", lbl, " | Stability Gain: ", damp*100, "%\n"))
  
  # Initialize State
  n_iter <- length(ms_risk$s_mu)
  mu_curr <- ms_risk$mu_trend[, 14] # Last observed trend (2023)
  mu_prev <- ms_risk$mu_trend[, 13]
  y_sim <- matrix(NA, nrow = n_iter, ncol = pred_years)
  
  for (t in 1:pred_years) {
    curr_year <- year_start + t - 1
    
    # Propagate Trend with Reduced Noise
    mu_next <- rnorm(n_iter, mean = 2 * mu_curr - mu_prev, sd = reduced_s_mu)
    
    # Inject Structural Shock in 2025
    if (curr_year == intervention_year) {
      mu_next <- mu_next + magnitude
    }
    
    # Generate Observations with Reduced Noise
    y_sim[, t] <- rnorm(n_iter, mean = mu_next, sd = reduced_s_Y)
    
    # Update State
    mu_prev <- mu_curr
    mu_curr <- mu_next
  }
  
  # Calculate Exceedance Probability for 2035
  y_2035 <- y_sim[, 12]
  p_exc <- mean(y_2035 > target_cap)
  
  # Store Results
  results_df <- bind_rows(results_df, data.frame(Value = y_2035, Scenario = lbl))
  summary_df <- bind_rows(summary_df, data.frame(Scenario = lbl, Level = lvl, P_exc = p_exc))
}

# 4. Print Summary Statistics
print(summary_df)

# 5. Plot Results (Ridge Plot)
results_df$Scenario <- factor(results_df$Scenario, levels = shock_labels)

p7_final <- ggplot(results_df, aes(x = Value, y = Scenario, fill = stat(x > target_cap))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha = 0.8) +
  geom_vline(xintercept = target_cap, linetype = "dashed", size = 1.2) +
  scale_fill_manual(values = c("TRUE" = "#de2d26", "FALSE" = "#2c7fb8"), 
                    labels = c("Compliant", "Non-Compliant"), name = "") +
  labs(
    title = "Regime Switch Analysis: From Drift to Anchor",
    x = "Log Consumption (2035)",
    y = ""
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

print(p7_final)
# ggsave("Figure7_RegimeSwitch.png", plot = p7_final, width = 8, height = 6)