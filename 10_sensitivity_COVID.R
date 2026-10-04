# =============================================================================
# 10_sensitivity_COVID.R
# COVID coding sensitivity (编码敏感性)：Primary / Pulse / Pulse-decay
# intput (输入)：Step2_primary_results.rds
# output (输出)：Step10_covid_sensitivity.rds
# =============================================================================
message("=== 10_sensitivity_COVID.R ===")
covid_schemes <- list(Primary= ifelse(CFG$years %in% c(2020, 2021, 2022), 1, 0),Pulse= ifelse(CFG$years == 2020, 1, 0),
                      Pulse_Decay = dplyr::case_when(CFG$years == 2020 ~ 1.00,CFG$years == 2021 ~ 0.50,
                                                     CFG$years == 2022 ~ 0.25,TRUE ~ 0))
classification_rows <- list()
diagnostic_rows     <- list()

for (scheme_idx in seq_along(covid_schemes)) {
  scheme_name  <- names(covid_schemes)[scheme_idx]
  covid_dummy  <- covid_schemes[[scheme_name]]
  cat("\nCOVID coding:", scheme_name, "\n")
  
  for (sector_idx in seq_along(SECTOR_CODES)) {
    sector <- SECTOR_CODES[sector_idx]
    pretty <- SECTOR_LABELS[[sector]]
    Y <- as.numeric(D[[sector]])
    
    full_data    <- make_full_data(Y, covid_dummy = covid_dummy)
    reduced_data <- make_reduced_data(Y, covid_dummy = covid_dummy)
    
    seed_base <- 200000 + scheme_idx * 10000 + sector_idx * 100
    lab <- paste(scheme_name, pretty)
    fit_full <- fit_retry(MODELS$Full_SORW, full_data, seed_base + 1, paste(lab, "Full"))
    fit_rsorw <- fit_retry(MODELS$Reduced_SORW, reduced_data, seed_base + 2, paste(lab, "Reduced SORW"))
    fit_rforw <- fit_retry(MODELS$Reduced_FORW, reduced_data, seed_base + 3, paste(lab, "Reduced FORW")) 
    
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- bind_rows(
      diagnose_stan_fit(fit_full,  "Full SORW",    pretty) %>% mutate(COVID_Coding = scheme_name),
      diagnose_stan_fit(fit_rsorw, "Reduced SORW", pretty) %>% mutate(COVID_Coding = scheme_name),
      diagnose_stan_fit(fit_rforw, "Reduced FORW", pretty) %>% mutate(COVID_Coding = scheme_name))
    
    loo_full  <- loo(extract_log_lik(fit_full,  "log_lik", merge_chains = FALSE), cores = 1)
    loo_rsorw <- loo(extract_log_lik(fit_rsorw, "log_lik", merge_chains = FALSE), cores = 1)
    loo_rforw <- loo(extract_log_lik(fit_rforw, "log_lik", merge_chains = FALSE), cores = 1)
    
    cmp_fr <- compare_pointwise(loo_full,  loo_rsorw, "Full SORW",    "Reduced SORW")
    cmp_rw <- compare_pointwise(loo_rforw, loo_rsorw, "Reduced FORW", "Reduced SORW")
    
    interpret_fr <- classify_evidence(cmp_fr$Delta_over_SE,"Full favored", "Reduced favored")
    interpret_rw <- classify_evidence(cmp_rw$Delta_over_SE,"FORW favored", "SORW favored")
    
    classification_rows[[length(classification_rows) + 1]] <- bind_rows(
      tibble(COVID_Coding = scheme_name, Sector = pretty,Comparison = "Full - Reduced", Delta_ELPD = cmp_fr$Delta_ELPD,
             SE_Delta = cmp_fr$SE_Delta, z = cmp_fr$Delta_over_SE,Interpretation = interpret_fr),
      tibble(COVID_Coding = scheme_name, Sector = pretty,Comparison = "FORW - SORW", Delta_ELPD = cmp_rw$Delta_ELPD,
             SE_Delta = cmp_rw$SE_Delta, z = cmp_rw$Delta_over_SE,Interpretation = interpret_rw))
    
    rm(fit_full, fit_rsorw, fit_rforw, loo_full, loo_rsorw, loo_rforw); gc()
  }
}

S10_model_comparison <- bind_rows(classification_rows)
S10_diagnostics      <- bind_rows(diagnostic_rows)

S10_classification_summary <- S10_model_comparison %>%select(COVID_Coding, Sector, Comparison, Interpretation) %>%
  pivot_wider(names_from = COVID_Coding, values_from = Interpretation)

covid_scheme_table <- tibble(Year = CFG$years,Primary = covid_schemes$Primary, Pulse   = covid_schemes$Pulse,
                             Pulse_Decay = covid_schemes$Pulse_Decay)

step10_covid <- list(covid_scheme_table = covid_scheme_table, model_comparison   = S10_model_comparison,
                     diagnostics= S10_diagnostics,classification= S10_classification_summary)
saveRDS(step10_covid, file.path(OUTPUT_DIR, "Step10_covid_sensitivity.rds"))

cat("10_sensitivity_COVID.R complete\n")