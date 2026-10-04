# =============================================================================
# 01_data_audit.R
# Data health checks. Does not modify the primary analysis. Outputs: DataAudit_*.csv
# (数据体检：不修改主分析，只输出 DataAudit_*.csv)
# Fixes the original execution-order bug: sector names are already defined in
# (修复原脚本执行顺序问题：sector 名称已在 00_setup.R 定义)
# =============================================================================
message("=== 01_data_audit.R ===")
# --- 1. Outcome ------------------------------------------------------
outcome_audit <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  y <- as.numeric(D[[sector]])
  log_jump <- diff(log(y))
  tibble(Sector= SECTOR_LABELS[[sector]], Min= min(y, na.rm = TRUE),Max= max(y, na.rm = TRUE),
         N_zero= sum(y <= 0, na.rm = TRUE),N_near_zero= sum(y < 1e-4, na.rm = TRUE),
         Min_positive= min(y[y > 0], na.rm = TRUE),Max_abs_log_change  = max(abs(log_jump), na.rm = TRUE),
         Year_max_abs_log_change = CFG$years[-1][which.max(abs(log_jump))],
         First_value= y[1],Last_value= y[length(y)])
})
#write_csv(outcome_audit, file.path(OUTPUT_DIR, "DataAudit_OutcomeIntegrity.csv"))

# --- 2. Year-specific --------------------------------------------------
outcome_long <-D%>%select(year,all_of(SECTOR_CODES))%>%pivot_longer(-year,names_to="Sector",values_to="Value")%>%
  mutate(Sector_Label = SECTOR_LABELS[Sector],Log_Value = log(Value),Log_Change = Log_Value - lag(Log_Value),
         Relative_Change = Value / lag(Value))

largest_jumps <- outcome_long %>%filter(!is.na(Log_Change)) %>%group_by(Sector_Label) %>%
  slice_max(abs(Log_Change), n = 3, with_ties = FALSE) %>%
  arrange(Sector_Label, desc(abs(Log_Change)))
#write_csv(largest_jumps, file.path(OUTPUT_DIR, "DataAudit_LargestOutcomeJumps.csv"))

# --- 3. LPG-specific audit (LPG 单独审计) --------------------------------------------------------
lpg_audit <- D %>% select(year, LiNGas) %>% mutate(Log_LPG = ifelse(LiNGas > 0, log(LiNGas), NA_real_),
                                                   Log_Change = Log_LPG - lag(Log_LPG),Ratio = LiNGas / lag(LiNGas))
#write_csv(lpg_audit, file.path(OUTPUT_DIR, "DataAudit_LPG.csv"))

# --- 4. Other-specific audit (Other 单独审计) ------------------------------------------------------
other_audit <- D %>% select(year, other) %>% mutate(Log_Other = ifelse(other > 0, log(other), NA_real_),
                                                    Log_Change = Log_Other - lag(Log_Other),Ratio = other / lag(other))
#write_csv(other_audit, file.path(OUTPUT_DIR, "DataAudit_Other.csv"))

# --- 5. Covariate scale (协变量尺度) ----------------------------------------------------------
X_raw    <- DIPAT %>% select(P, A, Tstr, Ttech)
X_log_df <- X_raw %>% mutate(across(everything(), log))
covariate_audit <- tibble(Variable = colnames(X_raw), Raw_Min= sapply(X_raw, min, na.rm = TRUE),
                          Raw_Max= sapply(X_raw, max, na.rm = TRUE),Raw_SD= sapply(X_raw, sd,  na.rm = TRUE),
                          Log_Min= sapply(X_log_df, min, na.rm = TRUE),Log_Max= sapply(X_log_df, max, na.rm = TRUE),
                          Log_SD= sapply(X_log_df, sd,  na.rm = TRUE),Log_Mean = sapply(X_log_df, mean, na.rm = TRUE))
#write_csv(covariate_audit, file.path(OUTPUT_DIR, "DataAudit_CovariateScale.csv"))

# --- 6. Correlation matrix (相关矩阵) ------------------------------------------------------------
cor_X <- cor(X_log_df, use = "pairwise.complete.obs", method = "pearson")
#write.csv(round(cor_X, 4), file.path(OUTPUT_DIR, "DataAudit_CovariateCorrelation.csv"))

# --- 7. VIF -----------------------------------------------------------------
calculate_vif <- function(df) {
  out <- numeric(ncol(df)); names(out) <- colnames(df)
  for (j in seq_along(out)) {
    y <- df[[j]]; others <- df[, -j, drop = FALSE]
    R2 <- summary(lm(y ~ ., data = others))$r.squared
    out[j] <- ifelse(is.finite(R2) && R2 < 1, 1 / (1 - R2), Inf)
  }
  out
}
vif_table <- tibble(Variable = colnames(X_log_df),VIF = as.numeric(calculate_vif(as.data.frame(X_log_df))))
#write_csv(vif_table, file.path(OUTPUT_DIR, "DataAudit_CovariateVIF.csv"))

# --- 8. Condition number ----------------------------------------------------
X_scaled <- scale(X_log_df)
cond_num <- kappa(X_scaled, exact = TRUE)
#writeLines(paste("Condition number:", cond_num),file.path(OUTPUT_DIR, "DataAudit_CovariateConditionNumber.txt"))

# --- 9. SD ratio flag -------------------------------------------------------
sd_ratio <- max(covariate_audit$Log_SD) / min(covariate_audit$Log_SD)
#writeLines(paste("SD ratio:", sd_ratio),file.path(OUTPUT_DIR, "DataAudit_CovariateSDRatio.txt"))

# ============================================================
# Additional diagnostic: correlation after removing linear time trend
# ============================================================
X_log_df <- DIPAT %>%select(P, A, Tstr, Ttech) %>%mutate(across(everything(), log))
time_index <- seq_len(nrow(X_log_df))
X_detrended <- as.data.frame(lapply(X_log_df, function(x) {resid(lm(x ~ time_index))}))
X_FULL <- as.matrix(X_detrended) 
colnames(X_detrended) <- colnames(X_log_df)
cor_X_detrended <- cor(X_detrended,use = "pairwise.complete.obs",method = "pearson")
print(round(cor_X_detrended, 3))
vif_detrended <- calculate_vif(X_detrended)
vif_detrended_table <- tibble(Variable = names(vif_detrended),VIF_detrended = as.numeric(vif_detrended))
cond_detrended <- kappa(scale(X_detrended), exact = TRUE)
cat("\nCondition number after linear detrending:",cond_detrended, "\n")


trend_R2 <- tibble(Variable = colnames(X_log_df),R2_YearTrend = sapply(
  X_log_df,function(x) summary(lm(x ~ time_index))$r.squared))
collinearity_audit <- covariate_audit %>%select(Variable, Log_SD) %>%left_join(vif_table, by = "Variable") %>%
  left_join(vif_detrended_table, by = "Variable") %>%left_join(trend_R2, by = "Variable")
collinearity_audit

cat("01_data_audit.R complete\n")