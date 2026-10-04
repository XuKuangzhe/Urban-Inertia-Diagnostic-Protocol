# =============================================================================
# 99_tables_figures.R (Publication-Ready Refactored Version)
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})
message("=== 99_tables_figures.R ===")
# --- 1. Load upstream results ------------------------------------------------
step2 <- readRDS(file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
step3 <- readRDS(file.path(OUTPUT_DIR, "Step3_identification.rds"))
step4 <- readRDS(file.path(OUTPUT_DIR, "Step4_forecast.rds"))
step6 <- if (file.exists(file.path(OUTPUT_DIR, "Step6_intervention.rds")))
  readRDS(file.path(OUTPUT_DIR, "Step6_intervention.rds")) else NULL
step7 <- if (file.exists(file.path(OUTPUT_DIR, "Step7_benchmark.rds")))
  readRDS(file.path(OUTPUT_DIR, "Step7_benchmark.rds")) else NULL
step8 <- if (file.exists(file.path(OUTPUT_DIR, "Step8_exact_LOO.rds")))
  readRDS(file.path(OUTPUT_DIR, "Step8_exact_LOO.rds")) else NULL
step10 <- if (file.exists(file.path(OUTPUT_DIR, "Step10_covid_sensitivity.rds")))
  readRDS(file.path(OUTPUT_DIR, "Step10_covid_sensitivity.rds")) else NULL

SECTORS_REV <- rev(SECTOR_ORDER_LABELS)
# =============================================================================
# Figure 1: Posterior Temporal-Persistence Diagnostic (R_pers)
# =============================================================================
p_Rpers <- ggplot(step3$persistence_table %>% mutate(Sector = factor(Sector, levels = SECTORS_REV)),aes(y = Sector)) +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "gray95", alpha = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = R_pers_Lower, xmax = R_pers_Upper), height = 0.2, linewidth = 0.8, color = "#1F4E79") +
  geom_point(aes(x = R_pers_Median), size = 3.5, color = "#1F4E79") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(title = "Posterior Temporal-Persistence Diagnostic",
       x = expression(R[pers] == sigma[mu]^2 / (sigma[mu]^2 + sigma[Y]^2)), y = NULL) + theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
#ggsave(file.path(FIG_DIR, "Figure1_Rpers.png"), p_Rpers, width = 7, height = 4.2, dpi = 300)

# =============================================================================
# Figure 2: Incremental Covariate Contribution (Exact-LOO Corrected)
# =============================================================================
fig2_data <- step8$comparison_full_reduced %>% mutate(Sector = factor(Sector, levels = SECTORS_REV),
                                                      CI_low  = Hybrid_Delta - 1.96 * Hybrid_SE,
                                                      CI_high = Hybrid_Delta + 1.96 * Hybrid_SE)

fig2 <- ggplot(fig2_data, aes(x = Hybrid_Delta, y = Sector)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18, linewidth = 0.8, color = "#C00000") +
  geom_point(size = 3.5, color = "#C00000") +
  #annotate("text", x = -2.8, y = 1.2, label = "Favors Reduced Model (Parsimony)", color = "gray30", size = 3.5, hjust = 0) +
  labs(title = "Incremental Predictive Contribution of Socioeconomic Covariates",
       x = expression(Delta * "ELPD"[LOO] ~ "(Full SORW - Reduced SORW)"), y = NULL) + theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(FIG_DIR, "Figure2_Incremental_Covariate.png"), fig2, width = 7.5, height = 4.2, dpi = 300)

# =============================================================================
# Figure 3: Predictive Evidence for Second-Order Temporal Dynamics
# =============================================================================
fig3_data <- pair_tab %>% filter(Basis == "hybrid") %>%
  mutate(Sector = factor(Sector, levels = SECTORS_REV), 
         Comparison = factor(Comparison, levels = c("FORW - SORW", "Drift - SORW", "Drift - FORW"),
                             labels = c("FORW - SORW", "DriftFORW - SORW", "DriftFORW - FORW")),
         CI_low  = Delta - 1.96 * SE, CI_high = Delta + 1.96 * SE)

fig3 <- ggplot(fig3_data, aes(x = Delta, y = Sector)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.8, color = "#2E75B6") +
  geom_point(size = 3, color = "#2E75B6") + facet_wrap(~ Comparison, ncol = 3, scales = "free_x") +
  labs(x = expression(Delta * "ELPD"), y = NULL) +theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Figure3_Temporal_Dynamics.png"), fig3, width = 10, height = 4.2, dpi = 300)


# =============================================================================
# Figure 4: Long-Horizon Exceedance & Policy Target Sensitivity
# =============================================================================
fig4a_data <- step4$analytic_2035_p_exc %>% 
  pivot_longer(c(P_exc_SORW, P_exc_FORW, P_exc_DRIFT), names_to = "Model", values_to = "P_exc") %>%
  mutate(Model = recode(Model, P_exc_SORW  = "SORW (2nd-order)", P_exc_FORW  = "FORW (1st-order)",
                        P_exc_DRIFT = "Drift-FORW (1st-order + drift)"),
         Sector = factor(Sector, levels = SECTOR_ORDER_LABELS))

fig4a <- ggplot(fig4a_data, aes(x = Sector, y = P_exc, fill = Model, shape = Model)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "black") +
  scale_shape_manual(values = c("SORW (2nd-order)" = 21, "FORW (1st-order)" = 24, "Drift-FORW (1st-order + drift)" = 22)) +
  scale_fill_manual(values = c("SORW (2nd-order)" = "#1F4E79", "FORW (1st-order)" = "#C00000", 
                               "Drift-FORW (1st-order + drift)" = "#2CA02C")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + coord_flip() +
  labs(title = "A: Exceedance Probability Under Baseline Target",
       x = NULL, y = expression(P(Y[2035] > 0.80 %*% Y[2023]))) + theme_bw(base_size = 11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))

fig4b_data <- step4$target_ratio_results %>% filter(Model == "SORW") %>%
  mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS))

fig4b <- ggplot(fig4b_data, aes(x = Target_Ratio, y = P_exc, color = Sector)) +
  geom_vline(xintercept = 0.80, linetype = "dotted", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0.50, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(0.5, 1.2, 0.1), labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "B: Sensitivity to Policy Target Strictness (SORW Momentum)",
       x = expression("Target Ratio Relative to 2023 Level (" * Y[2035] / Y[2023] * ")"),
       y = expression(P[exc]), color = "Sector") + theme_bw(base_size = 11) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 11))

fig4_combined <- fig4a / fig4b + plot_layout(heights = c(1, 1))
#ggsave(file.path(FIG_DIR, "Figure4_Target_Exceedance_and_Sensitivity.png"), fig4_combined, width = 8, height = 7.5, dpi = 300)
ggsave(file.path(FIG_DIR, "Figure4_Target_Exceedance.png"), fig4_combined, width = 8, height = 7, dpi = 300)

# =============================================================================
# Figure 5: Uncertainty Accumulation (5 部门整齐网格)
# =============================================================================
uncertainty_df <- step4$forecast_long %>%
  mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS),
         Model  = factor(Model, levels = c("SORW", "FORW"), labels = c("SORW (2nd-order)", "FORW (1st-order)")),
         Log_Interval_Width = Log_Upper_95 - Log_Lower_95)

p_uncertainty <- ggplot(uncertainty_df, aes(x = Year, y = Log_Interval_Width, color = Model, linetype = Model)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) + facet_wrap(~ Sector, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = seq(2024, 2035, by = 3)) +
  scale_color_manual(values = c("SORW (2nd-order)" = "#1F4E79", "FORW (1st-order)" = "#E77E23")) +
  labs(title = "Accumulation of Long-Horizon Predictive Uncertainty",
       x = "Forecast Year", y = expression(log(Q[97.5] / Q[2.5])),
       color = "State Specification", linetype = "State Specification") + theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "gray95"), plot.title = element_text(face = "bold"))
#ggsave(file.path(FIG_DIR, "Figure5_LongHorizon_Uncertainty.png"), p_uncertainty, width = 8.5, height = 5.5, dpi = 300)

# =============================================================================
# Figure 6: Intervention Dose-Response Surface (Shock + Damping 全貌)
# =============================================================================
if (!is.null(step6)) {
  fig6_data <- step6$results %>%
    mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS),
           Damping_Factor = factor(Damping, levels = c(0, 0.25, 0.50, 0.75),
                                   labels = c("Damping 0% (Shock only)", "Damping 25%", "Damping 50%", "Damping 75%")))
  
  p_intervention <- ggplot(fig6_data, aes(x = Reduction_pct, y = P_exc, color = Damping_Factor, group = Damping_Factor)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_line(linewidth = 0.85) + geom_point(size = 1.8) + facet_wrap(~ Sector, ncol = 3) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 50, 10)) + scale_color_brewer(palette = "Blues", direction = 1) +
    labs(title = "Conditional Reduction Response Surface",
         x = "Permanent Consumption Reduction Level from 2025 (%)",
         y = expression(P[exc]), color = "Volatility Mitigation") + theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.background = element_rect(fill = "gray95"), plot.title = element_text(face = "bold"))
  
  #ggsave(file.path(FIG_DIR, "Figure6_Intervention_DoseResponse.png"), p_intervention, width = 9, height = 6, dpi = 300)
}
ggsave(file.path(FIG_DIR, "Figure5_Intervention_Dose_Response.png"), p_intervention, width = 9, height = 6, dpi = 300)

# =============================================================================
# add S1 (SI): 4-Way Dynamic Benchmark LFO (预测评分对比)
# =============================================================================
if (!is.null(step7) && !is.null(step7$lfo_summary)) {
  p_lfo <- ggplot(step7$lfo_summary %>% mutate(Sector = factor(Sector, levels = SECTOR_ORDER_LABELS)),
                  aes(x = Model, y = Mean_LPD, fill = Model)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = Mean_LPD - SD_LPD/sqrt(N), ymax = Mean_LPD + SD_LPD/sqrt(N)),width = 0.2, linewidth = 0.6) +
    facet_wrap(~ Sector, ncol = 5, scales = "free_y") + scale_fill_brewer(palette = "Spectral") +
    labs(title = "Figure S1: Forward-Chaining Leave-Future-Out (LFO) Benchmark Scores",
         x = NULL, y = "Mean Log Predictive Density") + theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "gray95"), plot.title = element_text(face = "bold"))
  
  #ggsave(file.path(FIG_DIR, "FigureS1_LFO_Benchmark.png"), p_lfo, width = 10, height = 3.8, dpi = 300)
}

# -----------------------------------------------------------------------------
# Figure S2: COVID-coding sensitivity
# -----------------------------------------------------------------------------
if (!is.null(step10)) {
  plot_df <- step10$model_comparison %>%
    mutate(COVID_Coding = factor(COVID_Coding,levels = c("Primary", "Pulse", "Pulse_Decay"),
                                 labels = c("Primary", "Pulse", "Pulse-decay")),
           Comparison = factor(Comparison,levels = c("Full - Reduced", "FORW - SORW"),
                               labels = c("Full SORW - Reduced SORW", "FORW - SORW")),
           Sector = factor(Sector, levels = SECTOR_ORDER_LABELS),
           Lower = Delta_ELPD - 2 * SE_Delta, Upper = Delta_ELPD + 2 * SE_Delta)
  
  p_s2 <- ggplot(plot_df, aes(x = Sector, y = Delta_ELPD,color = COVID_Coding, group = COVID_Coding)) +
    geom_hline(yintercept = 0, linewidth = 0.6) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),position = position_dodge(width = 0.68),width = 0, linewidth = 0.65) +
    geom_point(position = position_dodge(width = 0.68), size = 2.8) +
    facet_wrap(~ Comparison, ncol = 1, scales = "fixed") +coord_flip() +
    scale_color_manual(values = c("Primary" = "#1F4E79","Pulse" = "#4F81BD","Pulse-decay" = "#A6CEE3")) +
    labs(x = NULL, y = expression(Delta*"ELPD"),color = "COVID coding") +theme_bw(base_size = 11)
}
ggsave(file.path(FIG_DIR, "FigureS2_COVID_sensitivity.png"), p_s2, width = 8.5, height = 7.5, dpi = 300)

# =============================================================================
# add S3 (SI)
# =============================================================================
p_coefs <- ggplot(step3$coefficient_table, aes(x = Median, y = Driver, color = Driver)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_errorbarh(aes(xmin = Lower_95, xmax = Upper_95), height = 0.25, linewidth = 0.7) +
  geom_point(size = 2.5) + facet_wrap(~ Sector, ncol = 5) +
  labs(title = "Figure S3: Posterior Elasticity of Socioeconomic Drivers (Full SORW)",
       x = "Regression Coefficient (Elasticity)", y = NULL) +
  scale_color_brewer(palette = "Dark2") + theme_bw(base_size = 11) +
  theme(legend.position = "none", strip.background = element_rect(fill = "gray95"), plot.title = element_text(face = "bold"))

#ggsave(file.path(FIG_DIR, "FigureS3_Socioeconomic_Coefficients.png"), p_coefs, width = 11, height = 3.5, dpi = 300)

# =============================================================================
# add S4 (SI)
# =============================================================================
fig4b_qc <- ggplot(step4$mc_vs_analytic %>%
                     pivot_longer(c(SORW_diff, FORW_diff), names_to = "Model", values_to = "Diff") %>%
                     mutate(Model = recode(Model, SORW_diff = "SORW", FORW_diff = "FORW")),
                   aes(x = Sector, y = Diff, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0) + scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Figure S4: Numerical Quality Control (Analytical P_exc vs Monte Carlo Simulation)",
       x = NULL, y = "Sampling Bias") + scale_fill_manual(values = c("SORW" = "#1F4E79", "FORW" = "#C00000")) +
  theme_bw(base_size = 11) + theme(legend.position = "top", plot.title = element_text(face = "bold"))

#ggsave(file.path(FIG_DIR, "FigureS4_MC_vs_Analytical_QC.png"), fig4b_qc, width = 7.5, height = 4, dpi = 300)

message("=== 99_tables_figures.R execution completed successfully. ===")

# -----------------------------------------------------------------------------
# Export tables (导出表格)
# -----------------------------------------------------------------------------
write_csv(step2$diagnostics,      file.path(OUTPUT_DIR, "Table_MCMC_Diagnostics.csv"))
write_csv(step2$Full_vs_Reduced,  file.path(OUTPUT_DIR, "Table_Full_vs_Reduced.csv"))
write_csv(step2$SORW_vs_FORW,     file.path(OUTPUT_DIR, "Table_SORW_vs_FORW.csv"))
write_csv(step3$persistence_table, file.path(OUTPUT_DIR, "Table_Rpers.csv"))
write_csv(step3$coefficient_table, file.path(OUTPUT_DIR, "Table_Coefficients.csv"))
write_csv(step4$analytic_2035_p_exc, file.path(OUTPUT_DIR, "Table_Pexc_2035.csv"))

library(knitr)
tabS_diag <- diagnostic_table %>% transmute(Carrier = Sector, Model = gsub("_", " ", Model),
                                            `Max Rhat` = round(Max_Rhat, 3), `Min ESS` = round(Min_n_eff),
                                            Divergences = Divergences, `Min BFMI` = round(Min_BFMI, 2))
writeLines(kable(tabS_diag, format = "latex", booktabs = TRUE, longtable = TRUE,
                 caption = "Sampling diagnostics of the 20 primary fits.", label = "si_sampling"),
           file.path(OUTPUT_DIR, "TabS_sampling_diagnostics.tex"))

cat("99_tables_figures.R complete\n")