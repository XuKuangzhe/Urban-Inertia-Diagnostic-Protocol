# ==============================================================================
# Script Name: Supplementary Material_Analysis.R
# Description: Generates Supplementary materials for validation and robustness checks.
#              Includes Parameter Recovery, Prior/Posterior Checks.
# Author: Xu Kuangzhe
# ==============================================================================

# ==============================================================================
# Part 0: Environment Check
# ==============================================================================
library(tidyverse)
library(rstan)
source("~/HDI.R")

# Ensure NullSSM model is compiled
if(!exists("NullSSM")) {
  stop("Error: 'NullSSM' is not found. Please compile 'SSM_Null_Dynamic.stan' first.")
}

# Ensure covid_vec exists (Default N=14, 2010-2023)
if(!exists("covid_vec")) {
  covid_vec <- rep(0, 14)
  names(covid_vec) <- 2010:2023
  covid_vec[names(covid_vec) %in% c("2020", "2021", "2022")] <- 1
}

# ==============================================================================
# Figure S2: BSTS Component Decomposition Plot (Null Model / Regime Model)
# ==============================================================================
target_energy_col <- "oil" #
energy_label <- "Oil Products"

cat(paste0("Generating Null Model Decomposition for: ", energy_label, "...\n"))

# --- 2. Null Model data ---
datalist_null_decomp <- list(
  T = 14, 
  Y = d[[target_energy_col]], 
  covid_dummy = covid_vec, 
  prior_alpha = 20, 
  prior_beta = 20
)

# --- 3. run Null Model  ---
fit_decomp_null <- sampling(NullSSM, data = datalist_null_decomp, iter = 4000, seed = 1234, refresh = 0)
post <- rstan::extract(fit_decomp_null)

# --- 4. ---
years <- 2010:2023

# A. Original data
obs_df <- data.frame(Year = years, Val = log(d[[target_energy_col]]), Type = "Observed")

# B. Inertial Trend: mu_trend
trend_mean <- colMeans(post$mu_trend)
trend_lower <- apply(post$mu_trend, 2, quantile, probs=0.025)
trend_upper <- apply(post$mu_trend, 2, quantile, probs=0.975)
df_trend <- data.frame(Year = years, Val = trend_mean, Lower = trend_lower, Upper = trend_upper)

# C. Exogenous Shock: beta_covid * D_t
beta_post <- post$beta_covid # Vector of posterior samples
shock_matrix <- matrix(NA, nrow = 4000, ncol = 14)
for(i in 1:4000){
  shock_matrix[i, ] <- beta_post[i] * covid_vec
}
shock_mean <- colMeans(shock_matrix)
shock_lower <- apply(shock_matrix, 2, quantile, probs=0.025)
shock_upper <- apply(shock_matrix, 2, quantile, probs=0.975)
df_shock <- data.frame(Year = years, Val = shock_mean, Lower = shock_lower, Upper = shock_upper)

# D. Total Fitted
total_mean <- colMeans(post$mu_total)

# --- 5. Plot (Null Model) ---
my_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

# sub plot 1: Data vs Fitted
p1 <- ggplot() +
  geom_point(data = obs_df, aes(x = Year, y = Val), size = 2, color = "black") +
  geom_line(aes(x = years, y = total_mean), color = "#08519c", size = 1) +
  labs(title = "A. Observed Data vs. Null Model Fit", y = "Log Consumption") +
  my_theme

# sub plot 2: Latent Inertia (Trend)
p2 <- ggplot(df_trend, aes(x = Year, y = Val)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#d95f0e", alpha = 0.2) +
  geom_line(color = "#d95f0e", size = 1) +
  labs(title = "B. Latent Inertial Trend (Shock-Free)", y = "Trend Level") +
  my_theme

# sub plot 3: Exogenous Shock (COVID-19)
p3 <- ggplot(df_shock, aes(x = Year, y = Val)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#7570b3", alpha = 0.2) +
  geom_line(color = "#7570b3", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "C. Exogenous Shock (COVID-19 Intervention)", y = "Impact Magnitude") +
  my_theme

# sub plot 4: Residuals
res_val <- obs_df$Val - total_mean
df_res <- data.frame(Year = years, Val = res_val)
p4 <- ggplot(df_res, aes(x = Year, y = Val)) +
  geom_col(fill = "gray40", width = 0.6) +
  geom_hline(yintercept = 0, color = "black") +
  labs(title = "D. Residuals (Unexplained Noise)", y = "Residual") +
  my_theme + theme(axis.title.x = element_text())

# --- 6. summary ---
final_plot_null <- (p1 / p2 / p3 / p4) + 
  plot_layout(heights = c(1, 1, 1, 0.8)) +
  plot_annotation(
    title = paste0("Null Model Decomposition: ", energy_label),
    subtitle = "Separating Latent Inertia from Exogenous Shocks."
  )

print(final_plot_null)
#plotpath<-paste0("~/Figure S2_Decomp_Null_", target_energy_col, ".png")
#ggsave(plotpath, plot = final_plot_null, width = 8, height = 10)


# ==============================================================================
# Figure S3: Prior vs. Posterior Comparison
# ==============================================================================

# 1. Run a representative model instance (e.g., Coal)
# ------------------------------------------------------------------------------
target_energy <- "coal" # Or use colnames(d)[2]
datalist_check <- list(T=14, Y=d[[target_energy]], covid_dummy=covid_vec, prior_alpha=20, prior_beta=20) 
fit_check <- sampling(NullSSM, data=datalist_check, iter=4000, seed=1234, refresh=0)

# 2. Extract Posterior Samples
# ------------------------------------------------------------------------------
ms <- rstan::extract(fit_check)
df_posterior <- bind_rows(
  data.frame(Value = ms$s_mu, Parameter = "s_mu (Trend Noise)", Type = "Posterior"),
  data.frame(Value = ms$s_Y,  Parameter = "s_Y (Obs Noise)",    Type = "Posterior")
)

# 3. Generate Theoretical Prior Distribution (Gamma 20, 20)
# ------------------------------------------------------------------------------
x_range <- seq(0, 2.5, length.out = 1000)
y_prior <- dgamma(x_range, shape = 20, rate = 20)

df_prior <- bind_rows(
  data.frame(Value = x_range, Density = y_prior, Parameter = "s_mu (Trend Noise)", Type = "Prior"),
  data.frame(Value = x_range, Density = y_prior, Parameter = "s_Y (Obs Noise)",    Type = "Prior")
)

# 4. Plot Comparison (Information Gain)
# ------------------------------------------------------------------------------
p_check <- ggplot() +
  # Posterior Density (Filled)
  geom_density(data = df_posterior, aes(x = Value, fill = Parameter), alpha = 0.6, color = NA) +
  # Prior Density (Dashed Line)
  geom_line(data = df_prior, aes(x = Value, y = Density, linetype = "Prior"), color = "black", linewidth = 0.8) +
  
  facet_wrap(~Parameter, scales = "free") +
  scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
  scale_linetype_manual(name = "", values = c("Prior" = "dashed")) +
  
  labs(title = "Prior vs. Posterior Comparison",
       subtitle = paste0("Information Gain Check (Sector: ", target_energy, ")"),
       x = "Parameter Value", y = "Density") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

print(p_check)
# ggsave("Figure S3_PriorPosterior.png", plot=p_check, width=8, height=5)

# ==============================================================================
# Figure S4: Sensitivity Analysis of Intervention Specification
# Description: Re-evaluates Model Comparison using a "Pulse-Decay" shock 
#              structure instead of a binary dummy variable.
# ==============================================================================

# 1. Define Pulse-Decay Intervention Vector
# ------------------------------------------------------------------------------
# Assumption: 
# 2020 = 100% shock magnitude
# 2021 = 50% (System adaptation)
# 2022 = 25% (Further adaptation)
# Others = 0
covid_vec_decay <- rep(0, 14)
names(covid_vec_decay) <- 2010:2023

covid_vec_decay["2020"] <- 1.0
covid_vec_decay["2021"] <- 0.5 
covid_vec_decay["2022"] <- 0.25

print("Pulse-Decay Vector defined:")
print(covid_vec_decay)

# 2. Run Model Comparison Loop
# ------------------------------------------------------------------------------
loo_summary_decay <- data.frame() 

# Ensure models are loaded
if(!exists("LSWaicM") | !exists("NullSSM")) {
  stop("Error: Stan models not found. Please compile them first.")
}

for(o in 2:8){
  energy_name <- colnames(d)[o] 
  print(paste("Processing Sensitivity (Pulse-Decay) for:", energy_name))
  
  # --- Data Pre-processing with Decay Vector ---
  # Replace binary covid_vec with covid_vec_decay in predictors
  X_log <- dipat[, 2:5] %>% mutate(covid_dummy = covid_vec_decay)
  X_log$P <- log(X_log$P)
  X_log$A <- log(X_log$A)
  X_log$Tstr <- log(X_log$Tstr)
  X_log$Ttech <- log(X_log$Ttech)
  
  # --- Prepare Stan Data Lists ---
  # Note: covid_dummy is updated to covid_vec_decay
  datalist_full <- list(T=14, K=5, X=X_log, Y=d[,o], prior_alpha=20, prior_beta=20)
  datalist_null <- list(T=14, Y=d[,o], covid_dummy=covid_vec_decay, prior_alpha=20, prior_beta=20) 
  
  # --- Sampling ---
  # Using iter=4000 to ensure convergence comparable to main analysis
  fit1 <- sampling(LSWaicM, data=datalist_full, iter=4000, seed=1234, refresh=0)
  fit2 <- sampling(NullSSM, data=datalist_null, iter=4000, seed=1234, refresh=0)
  
  # --- LOO Calculation ---
  log_lik_full <- extract_log_lik(fit1, merge_chains = FALSE)
  log_lik_null <- extract_log_lik(fit2, merge_chains = FALSE)
  
  loo_full <- loo(log_lik_full)
  loo_null <- loo(log_lik_null)
  
  # --- Comparison ---
  comp <- loo_compare(list(Full_Model = loo_full, Null_Model = loo_null))
  
  # --- Store Results ---
  comp_df <- as.data.frame(comp) %>%
    rownames_to_column(var = "Model_Name") %>%
    mutate(Energy = energy_name)
  
  loo_summary_decay <- rbind(loo_summary_decay, comp_df)
}

# 3. Plotting
# ------------------------------------------------------------------------------
plot_data_decay <- loo_summary_decay %>%
  filter(Model_Name == "Full_Model") %>%
  mutate(
    Energy_Label = recode(Energy, !!!energy_labels), # Use the standard mapping list
    Energy_Label = reorder(Energy_Label, elpd_diff),
    lower = elpd_diff - 2 * se_diff,
    upper = elpd_diff + 2 * se_diff,
    is_significant = upper < 0 
  )

p_decay <- ggplot(plot_data_decay, aes(x = elpd_diff, y = Energy_Label)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_col(aes(fill = is_significant), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.6) +
  
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#99C9E7"), guide = "none") +
  
  labs(
    title = "Sensitivity Analysis of Model Selection",
    subtitle = "Robustness check using Pulse-Decay Intervention (1.0, 0.5, 0.25)",
    x = expression(paste(Delta, "ELPD (Full Model - Null Model)")),
    y = ""
  ) +
  
  annotate("text", x = -0.5, y = 0.5, label = "Favors Null Model\n(Inertia)", 
           hjust = 1, size = 3.5, fontface = "italic", color = "gray30") +
  
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major.y = element_blank())

print(p_decay)

# Save the plot
# ggsave("Figure S4_Sensitivity_Intervention.png", plot=p_decay, width=9, height=5)


# ==============================================================================
# Table S1: Coefficient Analysis
# ==============================================================================
# Objective: Extract regression coefficients from the Full Model (STIRPAT) 
# to demonstrate that the posterior intervals for all drivers include zero.
# This generates the data required for Appendix Table B.1.

coef_summary_df <- data.frame()

# Define STIRPAT variable names (corresponding to the first 4 columns of X_log)
# Note: In datalist_full, K=5 and the 5th column is the Covid dummy. 
# We focus primarily on the first 4 socio-economic drivers.
var_names <- c("Population (P)", "Affluence (A/GDP)", "Structure (Str)", "Technology (Tech)")

# Loop through each energy sector
for(o in 2:8){
  energy_name <- colnames(d)[o]
  print(paste("Extracting Coefficients for:", energy_name))
  
  # --- 1. Data Preparation (Consistent with Part 1) ---
  X_log <- dipat[, 2:5] %>% mutate(covid_dummy = covid_vec)
  X_log$P <- log(X_log$P)
  X_log$A <- log(X_log$A)
  X_log$Tstr <- log(X_log$Tstr)
  X_log$Ttech <- log(X_log$Ttech)
  
  datalist_full <- list(T=14, K=5, X=X_log, Y=d[,o], prior_alpha=20, prior_beta=20)
  
  # --- 2. Run Full Model ---
  # Ideally, if fit1 was saved in Part 1, reuse it here to save time.
  # For clarity in this standalone script, we re-run sampling (approx. 1-2 mins/sector).
  fit1 <- sampling(LSWaicM, data=datalist_full, iter=4000, seed=1234, refresh=0)
  
  # --- 3. Extract Beta Coefficients ---
  # Assuming the parameter name for regression coefficients in the Stan model is "beta".
  # Returns a matrix of [iterations, K].
  post_beta <- rstan::extract(fit1)$beta 
  
  # --- 4. Calculate Statistics ---
  for(k in 1:4){ # Iterate through the first 4 STIRPAT drivers (skipping Covid)
    beta_samples <- post_beta[, k]
    
    # Calculate Mean and HDI (Highest Density Interval)
    # Ensure the summary_MCMC function is defined in your environment
    b_mean <- summary_MCMC(beta_samples)["mean"]
    b_hdi_h <- summary_MCMC(beta_samples)["HDI_high"]
    b_hdi_l <- summary_MCMC(beta_samples)["HDI_low"]
    
    # Calculate Probability of Direction (pd)
    # pd is the probability that the effect goes in the direction of the sign.
    # Values close to 0.5 indicate randomness; 1.0 indicates certainty.
    # We also check if zero is contained within the HDI (Zero_in_HDI).
    prob_pos <- mean(beta_samples > 0)
    prob_direction <- max(prob_pos, 1 - prob_pos) 
    
    # Compile results
    row <- data.frame(
      Sector = energy_name,
      Driver = var_names[k],
      Mean   = round(b_mean, 3),
      Lower_95 = round(b_hdi_l, 3),
      Upper_95 = round(b_hdi_h, 3),
      Zero_in_HDI = (b_hdi_l < 0 & b_hdi_h > 0), # TRUE indicates structural insignificance
      Prob_Direction = round(prob_direction, 2)
    )
    coef_summary_df <- rbind(coef_summary_df, row)
  }
}

# --- 5. Export Table ---
# Recode sector names for readability using predefined labels
coef_summary_df <- coef_summary_df %>%
  mutate(Sector = recode(Sector, !!!energy_labels))

print("=== Appendix Table: STIRPAT Coefficients Analysis ===")
print(coef_summary_df)

# Save to CSV
write_csv(coef_summary_df, "Table_Appendix_Coefficients.csv")

# --- 6. (Optional) Visualization of Coefficients ---
# This plot can be included in the appendix to visually support the argument.
p_coef <- ggplot(coef_summary_df, aes(x = Mean, y = Driver)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = Lower_95, xmax = Upper_95), height = 0.2) +
  geom_point(size = 3, color = "#2c7fb8") +
  facet_wrap(~Sector, scales = "free_x") +
  labs(
    title = "Why the STIRPAT Model Failed: Coefficient Analysis",
    subtitle = "Posterior Means & 95% HDI for Socio-economic Drivers",
    x = "Elasticity Coefficient", y = ""
  ) +
  theme_bw(base_size = 16)

print(p_coef)
# ggsave("Figure_Appendix_Coefficients.png", plot=p_coef, width=10, height=8)