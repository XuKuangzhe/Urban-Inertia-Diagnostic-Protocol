# =============================================================================
# 02_primary_models.R
# PRIMARY: 3 models × 7 sectors
#   fitting + MCMC diagnostics + PSIS-LOO + model comparisons
# Output:  Step2_primary_results.rds
# =============================================================================
message("=== 02_primary_models.R ===")
# =============================================================================
# 0. Storage objects
# =============================================================================
fits <- list()
loo_objects <- list()
diagnostic_results <- list()
comparison_results <- list()
fit_metadata <- list()
# =============================================================================
# 1. Model-fitting helper with automatic retry
# =============================================================================
fit_model_with_retry <- function(model, data, seed, model_name, sector_name,iter = CFG$iter,
                                 adapt_delta_grid = c(CFG$adapt_delta, 0.995,0.999),
                                 max_treedepth = CFG$max_treedepth) {
  adapt_delta_grid <- unique(adapt_delta_grid)
  retry_log <- list()
  
  for (i in seq_along(adapt_delta_grid)) {
    this_delta <- adapt_delta_grid[i]
    cat("\n>>> Fitting:", sector_name, "-", model_name,
        "| adapt_delta =", this_delta,
        "| attempt =", i, "/", length(adapt_delta_grid), "\n")
    fit <- tryCatch(
      sampling(model, data = data, iter = iter, seed = seed,
               control = list(adapt_delta = this_delta, max_treedepth = max_treedepth),refresh = 0),
      error = function(e) structure(list(error_message = conditionMessage(e)),class = "stan_fit_error")
    )
    # --- Sampling itself failed ---
    if (inherits(fit, "stan_fit_error")) {
      retry_log[[i]] <- tibble(Attempt = i, Adapt_Delta = this_delta,
                               Sampling_Error = TRUE, Error_Message = fit$error_message)
      cat("Sampling error:", fit$error_message, "\n")
      next
    }
    # --- Diagnose fitted model ---
    diag <- diagnose_stan_fit(fit, model_name, sector_name)
    retry_log[[i]] <- tibble(Attempt = i, Adapt_Delta = this_delta, Sampling_Error = FALSE,
                             Error_Message = NA_character_,Max_Rhat = diag$Max_Rhat, 
                             Min_n_eff = diag$Min_n_eff,Divergences = diag$Divergences, 
                             Max_Treedepth_Hits = diag$Max_Treedepth_Hits,Min_BFMI = diag$Min_BFMI,
                             Rhat_OK = diag$Rhat_OK, Divergence_OK = diag$Divergence_OK,
                             Treedepth_OK = diag$Treedepth_OK, BFMI_OK = diag$BFMI_OK)
    # --- Model passes all diagnostics ---
    all_ok <- diag$Rhat_OK && diag$Divergence_OK && diag$Treedepth_OK && diag$BFMI_OK
    if (all_ok) {
      cat("PASS:", sector_name, "-", model_name, "| adapt_delta =", this_delta, "\n")
      return(list(
        fit = fit, diagnostics = diag,
        metadata = list(seed = seed, iter = iter, adapt_delta_used = this_delta,
                        max_treedepth = max_treedepth, retry_number = i - 1,retry_log = bind_rows(retry_log))))
    }
    # --- Diagnostic failure: prepare next retry ---
    cat("Diagnostic failure:", sector_name, "-", model_name,
        "| divergences =", diag$Divergences,
        "| max Rhat =", round(diag$Max_Rhat, 5), "\n")
  }
  
  retry_table <- bind_rows(retry_log)
  print(retry_table)
  
  stop(paste0(
    "\nPRIMARY MODEL FAILURE\n",
    "Sector: ", sector_name, "\n",
    "Model: ", model_name, "\n",
    "No acceptable posterior sample was obtained after ",
    length(adapt_delta_grid), " attempts.\n",
    "Primary analysis has been stopped rather than calculating LOO from ",
    "an inadequately diagnosed fit."))
}

# =============================================================================
# 2. Main sector loop
# =============================================================================
for (sector in SECTOR_CODES) {
  pretty <- SECTOR_LABELS[[sector]]
  cat("\n============================================================\n")
  cat("Sector:", pretty, "\n")
  cat("============================================================\n")
  Y <- as.numeric(D[[sector]])
  stopifnot(length(Y) == length(CFG$years), all(is.finite(Y)), all(Y > 0))
  full_data <- make_full_data(Y)
  reduced_data <- make_reduced_data(Y)
  sector_index <- match(sector, SECTOR_CODES)
  seed_s <- CFG$seed + 100 * sector_index
  seeds <- c(Full_SORW = seed_s + 1, Reduced_SORW = seed_s + 2, Reduced_FORW = seed_s + 3, Reduced_DriftFORW = seed_s + 4)
  fit_full   <- fit_model_with_retry(MODELS$Full_SORW,         full_data,    seeds["Full_SORW"],         "Full SORW",         pretty)
  fit_rsorw  <- fit_model_with_retry(MODELS$Reduced_SORW,      reduced_data, seeds["Reduced_SORW"],      "Reduced SORW",      pretty)
  fit_rforw  <- fit_model_with_retry(MODELS$Reduced_FORW,      reduced_data, seeds["Reduced_FORW"],      "Reduced FORW",      pretty)
  fit_rdrift <- fit_model_with_retry(MODELS$Reduced_DriftFORW, reduced_data, seeds["Reduced_DriftFORW"], "Reduced DriftFORW", pretty)
  
  fits[[sector]] <- list(Full_SORW = fit_full$fit, Reduced_SORW = fit_rsorw$fit,
                         Reduced_FORW = fit_rforw$fit, Reduced_DriftFORW = fit_rdrift$fit)
  sector_diag <- bind_rows(fit_full$diagnostics, fit_rsorw$diagnostics,fit_rforw$diagnostics, fit_rdrift$diagnostics)
  diagnostic_results[[sector]] <- sector_diag
  fit_metadata[[sector]] <- list(Full_SORW = fit_full$metadata, Reduced_SORW = fit_rsorw$metadata,
                                 Reduced_FORW = fit_rforw$metadata, Reduced_DriftFORW = fit_rdrift$metadata)
  
  ll_full   <- extract_log_lik(fit_full$fit,   "log_lik", merge_chains = FALSE)
  ll_rsorw  <- extract_log_lik(fit_rsorw$fit,  "log_lik", merge_chains = FALSE)
  ll_rforw  <- extract_log_lik(fit_rforw$fit,  "log_lik", merge_chains = FALSE)
  ll_rdrift <- extract_log_lik(fit_rdrift$fit, "log_lik", merge_chains = FALSE)
  loo_full   <- loo(ll_full,   cores = 1)
  loo_rsorw  <- loo(ll_rsorw,  cores = 1)
  loo_rforw  <- loo(ll_rforw,  cores = 1)
  loo_rdrift <- loo(ll_rdrift, cores = 1)
  
  loo_objects[[sector]] <- list(Full_SORW = loo_full, Reduced_SORW = loo_rsorw,
                                Reduced_FORW = loo_rforw, Reduced_DriftFORW = loo_rdrift)
  
  comp_fr <- compare_pointwise(loo_full,   loo_rsorw, "Full_SORW",         "Reduced_SORW")
  comp_rw <- compare_pointwise(loo_rforw,  loo_rsorw, "Reduced_FORW",      "Reduced_SORW")
  comp_dr <- compare_pointwise(loo_rdrift, loo_rsorw, "Reduced_DriftFORW", "Reduced_SORW")
  
  comparison_results[[sector]] <- list(Full_vs_Reduced = comp_fr, FORW_vs_SORW = comp_rw, Drift_vs_SORW = comp_dr)
  
  checkpoint <- list(fits = fits, loo = loo_objects, diagnostics = diagnostic_results,
                     comparisons = comparison_results, fit_metadata = fit_metadata)
  saveRDS(checkpoint, file.path(OUTPUT_DIR, "Step2_primary_checkpoint.rds"))
  gc(verbose = FALSE)
  cat("\nSector completed:", pretty, "\n")
}
# =============================================================================
# 3. Summary diagnostic table
# =============================================================================
diagnostic_table <- bind_rows(diagnostic_results) %>% select(Sector, Model, everything())
# =============================================================================
# 4. Full vs Reduced summary
# =============================================================================
full_reduced_table <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  loo_full <- loo_objects[[sector]]$Full_SORW
  loo_red  <- loo_objects[[sector]]$Reduced_SORW
  cmp <- comparison_results[[sector]]$Full_vs_Reduced
  
  tibble(Sector = SECTOR_LABELS[[sector]],ELPD_Full = loo_full$estimates["elpd_loo", "Estimate"],
         ELPD_Reduced = loo_red$estimates["elpd_loo", "Estimate"],
         Delta_ELPD_Full_minus_Reduced = cmp$Delta_ELPD, SE_Delta = cmp$SE_Delta,
         Delta_over_SE = cmp$Delta_over_SE,
         Evidence = classify_evidence(cmp$Delta_over_SE,pos_label = "Full favored",neg_label = "Reduced favored"))
})

# =============================================================================
# 5. FORW vs SORW summary
# =============================================================================
rw_table <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  loo_sorw <- loo_objects[[sector]]$Reduced_SORW
  loo_forw <- loo_objects[[sector]]$Reduced_FORW
  cmp <- comparison_results[[sector]]$FORW_vs_SORW
  
  tibble(Sector = SECTOR_LABELS[[sector]],ELPD_SORW = loo_sorw$estimates["elpd_loo", "Estimate"],
         ELPD_FORW = loo_forw$estimates["elpd_loo", "Estimate"],Delta_ELPD_FORW_minus_SORW = cmp$Delta_ELPD,
         SE_Delta = cmp$SE_Delta,Delta_over_SE = cmp$Delta_over_SE,
         Evidence = classify_evidence(cmp$Delta_over_SE,pos_label = "FORW favored",neg_label = "SORW favored"))
})

drift_table <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  cmp <- comparison_results[[sector]]$Drift_vs_SORW
  tibble(Sector = SECTOR_LABELS[[sector]],
         ELPD_SORW = loo_objects[[sector]]$Reduced_SORW$estimates["elpd_loo", "Estimate"],
         ELPD_Drift = loo_objects[[sector]]$Reduced_DriftFORW$estimates["elpd_loo", "Estimate"],
         Delta_ELPD_Drift_minus_SORW = cmp$Delta_ELPD, SE_Delta = cmp$SE_Delta,
         Delta_over_SE = cmp$Delta_over_SE,
         Evidence = classify_evidence(cmp$Delta_over_SE, "Drift-FORW favored", "SORW favored"))
})

write_csv(drift_table, file.path(OUTPUT_DIR, "Step2_Drift_vs_SORW_LOO.csv"))
# =============================================================================
# 6. Pareto-k audit
# =============================================================================
pareto_audit <- purrr::map_dfr(SECTOR_CODES, function(sector) {
  purrr::map_dfr(names(loo_objects[[sector]]), function(model_name) {
    loo_obj <- loo_objects[[sector]][[model_name]]
    k <- as.numeric(pareto_k_values(loo_obj))
    tibble(Sector = sector, Pretty = SECTOR_LABELS[[sector]], Model = model_name,
           Year = CFG$years, Index = seq_along(CFG$years), Pareto_k = k,
           Flag = case_when(k > 1.0 ~ ">1.0 HIGH PRIORITY",k > 0.7 ~ ">0.7 CHECK",TRUE    ~ "OK"))
  })
}) %>% arrange(desc(Pareto_k))

# =============================================================================
# 7. Pareto-k summary
# =============================================================================
pareto_summary <- pareto_audit %>%group_by(Sector, Pretty, Model) %>%
  summarise(Max_Pareto_k = max(Pareto_k),Num_k_gt_0_7 = sum(Pareto_k > 0.7),
            Num_k_gt_1_0 = sum(Pareto_k > 1.0),.groups = "drop")

# =============================================================================
# 8. Final primary result object
# =============================================================================
step2_results <- list(metadata = list(analysis = "Primary Bayesian structural time-series analysis",
                                      years = CFG$years, sectors = SECTOR_CODES,
                                      models = c("Full_SORW", "Reduced_SORW", "Reduced_FORW", "Reduced_DriftFORW"),
                                      prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta,
                                      lasso_alpha = CFG$lasso_alpha, lasso_beta = CFG$lasso_beta,
                                      adapt_delta_initial = CFG$adapt_delta,
                                      adapt_delta_retry = c(CFG$adapt_delta, 0.995, 0.999),
                                      max_treedepth = CFG$max_treedepth, n_iter = CFG$iter,
                                      seed_scheme = "sector-specific seed = CFG$seed + 100 * sector index; "),
                      fits = fits, loo = loo_objects, diagnostics = diagnostic_table,
                      Full_vs_Reduced = full_reduced_table, SORW_vs_FORW = rw_table, Drift_vs_SORW = drift_table,
                      pareto_audit = pareto_audit, pareto_summary = pareto_summary, fit_metadata = fit_metadata)

# =============================================================================
# 9. Save final primary objects
# =============================================================================
saveRDS(step2_results, file.path(OUTPUT_DIR, "Step2_primary_results.rds"))
write_csv(diagnostic_table,  file.path(OUTPUT_DIR, "Step2_MCMC_Diagnostics.csv"))
write_csv(full_reduced_table, file.path(OUTPUT_DIR, "Step2_Full_vs_Reduced_LOO.csv"))
write_csv(rw_table,          file.path(OUTPUT_DIR, "Step2_SORW_vs_FORW_LOO.csv"))
write_csv(pareto_audit,      file.path(OUTPUT_DIR, "Step2_Pareto_k_Audit.csv"))
write_csv(pareto_summary,    file.path(OUTPUT_DIR, "Step2_Pareto_k_Summary.csv"))

# =============================================================================
# 10. Final status check
# =============================================================================
stopifnot(length(fits) == length(SECTOR_CODES), length(loo_objects) == length(SECTOR_CODES),
          length(diagnostic_results) == length(SECTOR_CODES),length(comparison_results) == length(SECTOR_CODES))

cat("\n============================================================\n")
cat("02_primary_models.R COMPLETE\n")
cat("Sectors:", length(fits), "/", length(SECTOR_CODES), "\n")
cat("LOO objects:", length(loo_objects), "/", length(SECTOR_CODES), "\n")
cat("============================================================\n")