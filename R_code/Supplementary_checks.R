# Supplementary checks
# =============================================================================
# supplementary_checks.R
# 在跑完 summaryAllStepCode.R 的【同一个 R 会话】里运行
# 需要的对象: step2_results, exact_audit_results, step2_data_v3, MODELS, CFG, D,
#             COVID_PRIMARY, SECTOR_CODES, SECTOR_LABELS, OUTPUT_DIR,
#             fit_prefix_model, log_mean_exp, step4_forecast
# 每个 section 独立，可单独运行。结果同时写入 OUTPUT_DIR/supp/ 下的 csv。
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(loo); library(rstan) })
stopifnot(exists("step2_results"), exists("CFG"), exists("MODELS"))
OUT <- file.path(OUTPUT_DIR, "supp"); dir.create(OUT, showWarnings = FALSE)

# =============================================================================
# A. 收敛与重试记录（对应问题 3）
# =============================================================================
## A1. Step 2 主模型的 adapt_delta 重试记录
retry_step2 <- purrr::imap_dfr(step2_results$fit_metadata, function(lst, sector) {
  purrr::imap_dfr(lst, function(m, model) {
    tibble(Sector = sector, Model = model,
           Retry_number = m$retry_number, Adapt_delta_used = m$adapt_delta_used)
  })
})
cat("\n[A1] Step 2: 各 adapt_delta / 重试次数的模型数\n")
print(dplyr::count(retry_step2, Adapt_delta_used, Retry_number))
cat("\n[A1] Step 2: 发生过重试的模型\n")
print(filter(retry_step2, Retry_number > 0), n = Inf)
#write_csv(retry_step2, file.path(OUT, "A1_retry_step2.csv"))

## A2. Step 9 / 9B / 10 的重试记录（fit_retry 把日志写进全局 RETRY_LOG）
cat("\n[A2] Step 9/9B/10 重试记录\n")
if (exists("RETRY_LOG") && length(RETRY_LOG) > 0) {
  retry_s9 <- bind_rows(RETRY_LOG)
  print(dplyr::count(retry_s9, Attempt_Used, Passed))
  print(filter(retry_s9, !Passed | Attempt_Used > 1), n = Inf)   # Passed=FALSE 表示三次都没过，返回的是"最不差"的拟合
  #write_csv(retry_s9, file.path(OUT, "A2_retry_step9_10.csv"))
} else message("RETRY_LOG 不在当前会话；请改发下面三个对象里 Divergence_OK/Rhat_OK/BFMI_OK 为 FALSE 的行")
for (nm in c("step9_prior", "step9B_lasso", "step10_covid")) {
  if (exists(nm)) {
    d <- get(nm)$diagnostics
    cat("\n[A2] ", nm, ": 未通过的拟合\n", sep = "")
    print(as.data.frame(d %>% filter(!Rhat_OK | !Divergence_OK | !BFMI_OK)))
  }
}

## A3. exact-LOO refit 数量与收敛
stopifnot(exists("exact_audit_results"))
cat("\n[A3] 期望 refit 数 = k>0.7 的点数 =", sum(step2_results$pareto_summary$Num_k_gt_0_7),
    "; 实际 refit 数 =", nrow(exact_audit_results), "\n")
exact_summary <- exact_audit_results %>% group_by(Model) %>%
  summarise(N_refits = n(), Max_Rhat = max(Max_Rhat), Min_n_eff = min(Min_n_eff),
            Max_Divergences = max(Divergences), Max_Div_Rate = max(Div_Rate),
            Max_Treedepth_hits = max(Treedepth_hits), N_not_OK = sum(!Convergence_OK), .groups = "drop")
print(as.data.frame(exact_summary))
print(as.data.frame(dplyr::count(exact_audit_results, Sector, Model)))
print(as.data.frame(filter(exact_audit_results, !Convergence_OK)))
#write_csv(exact_audit_results, file.path(OUT, "A3_exact_refit_log.csv"))

## A4. 03B: R_pers 先验 vs 后验，并检查 sigma 的后验是否落在先验 95% 区间之外
cat("\n[A4] R_pers prior vs posterior (03B)\n")
print(as.data.frame(Rpers_prior_posterior))
cat("prior R_pers 2.5/50/97.5% =", round(prior_q, 3), "; prior P(R_pers>0.5) =",round(mean(prior_Rpers > 0.5), 3), "\n")
pr_mu <- qlnorm(c(.025, .5, .975), -2.0, 0.7); pr_Y <- qlnorm(c(.025, .5, .975), -2.3, 0.6)
sigma_pp <- purrr::map_dfr(SECTOR_CODES, function(s) {
  p <- rstan::extract(step2_results$fits[[s]]$Reduced_SORW, pars = c("s_mu", "s_Y"))
  tibble(Sector = SECTOR_LABELS[[s]],
         s_mu_post_med = median(p$s_mu), s_mu_post_lo = quantile(p$s_mu, .025), s_mu_post_hi = quantile(p$s_mu, .975),
         s_Y_post_med = median(p$s_Y),  s_Y_post_lo = quantile(p$s_Y, .025),  s_Y_post_hi = quantile(p$s_Y, .975),
         s_mu_med_below_prior_q025 = median(p$s_mu) < pr_mu[1],
         s_Y_med_below_prior_q025  = median(p$s_Y)  < pr_Y[1])
})
cat("prior sigma_mu 2.5/50/97.5% =", round(pr_mu, 4), "; prior sigma_Y =", round(pr_Y, 4), "\n")
print(as.data.frame(sigma_pp))
#write_csv(sigma_pp, file.path(OUT, "A4_sigma_prior_posterior.csv"))

# =============================================================================
# B. v2 与 v3 (use_obs 全为 1) 是否等价（Full SORW）——先验 Jacobian 检查
#    v3 里 `beta ~ double_exponential(0, S_beta)` 作用在变换后的参数上、没有 Jacobian，
#    等价于多乘了 1/S_beta，即 S_beta 先验由 Gamma(30,30) 变为 Gamma(29,30)（均值 1 -> 0.967）。
#    修复：v3 里把该行改成  beta_raw ~ double_exponential(0, 1);  （与 v2 一致）
# =============================================================================
chk <- "electric"
fit_v3 <- sampling(MODELS$Full_SORW_v3, data = step2_data_v3[[chk]]$Full_SORW,iter = CFG$iter, seed = 2026,
                   control = list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth), refresh = 0)
fit_v2 <- step2_results$fits[[chk]]$Full_SORW
tab_fit <- function(f, lab) {
  s <- summary(f, pars = c("S_beta", "s_mu", "s_Y", "beta_covid"))$summary
  tibble(Version = lab, Par = rownames(s), Mean = s[, "mean"], SD = s[, "sd"], n_eff = s[, "n_eff"])
}
cat("\n[B] v2 vs v3 (", chk, ")  期望: 仅 S_beta 的均值低约 0.03 (先验 1 -> 0.967)\n", sep = "")
print(as.data.frame(bind_rows(tab_fit(fit_v2, "v2"), tab_fit(fit_v3, "v3"))))

# =============================================================================
# C. 逐点结果导出 + Drift vs FORW (hybrid) + stacking 权重（对应问题 4）
# =============================================================================
models4 <- c("Full_SORW", "Reduced_SORW", "Reduced_FORW", "Reduced_DriftFORW")
pw_long <- purrr::map_dfr(SECTOR_CODES, function(s) {
  purrr::map_dfr(models4, function(m) {
    pw <- as.data.frame(step2_results$loo[[s]][[m]]$pointwise)
    tibble(Sector = s, Model = m, Index = seq_len(nrow(pw)), Year = CFG$years, elpd_psis = pw$elpd_loo)
  })
}) %>%
  left_join(step2_results$pareto_audit %>% select(Sector, Model, Index, Pareto_k), by = c("Sector", "Model", "Index")) %>%
  left_join(exact_audit_results %>% select(Sector, Model, Index, elpd_exact = Exact_ELPD), by = c("Sector", "Model", "Index")) %>%
  mutate(elpd_hybrid = coalesce(elpd_exact, elpd_psis))
#write_csv(pw_long, file.path(OUT, "C_pointwise_elpd_hybrid.csv"))  
hi <- pw_long %>% filter(Pareto_k > 0.7) %>% mutate(delta = elpd_exact - elpd_psis)
c(neg = sum(hi$delta < 0), median = median(hi$delta), mean = mean(hi$delta),
  gt05 = sum(abs(hi$delta) > 0.5), gt1 = sum(abs(hi$delta) > 1), maxneg = min(hi$delta))
hi %>% group_by(Sector) %>% summarise(mean = mean(delta))
hi %>% dplyr::count(Year)

cat("\n[C] 已写出 C_pointwise_elpd_hybrid.csv, 行数 =", nrow(pw_long), "(应为 280)\n")

wide_h <- function(s, col) pw_long %>% filter(Sector == s) %>% select(Index, Model, all_of(col)) %>%
  pivot_wider(names_from = Model, values_from = all_of(col)) %>% arrange(Index)
cmp <- function(a, b) { d <- a - b; n <- length(d); tibble(Delta = sum(d), SE = sqrt(n * var(d)), z = sum(d) / sqrt(n * var(d))) }

pair_tab <- purrr::map_dfr(SECTOR_CODES, function(s) {
  purrr::map_dfr(c("elpd_psis", "elpd_hybrid"), function(col) {
    w <- wide_h(s, col)
    bind_rows(
      cmp(w$Reduced_DriftFORW, w$Reduced_FORW) %>% mutate(Comparison = "Drift - FORW"),
      cmp(w$Reduced_DriftFORW, w$Reduced_SORW) %>% mutate(Comparison = "Drift - SORW"),
      cmp(w$Reduced_FORW,      w$Reduced_SORW) %>% mutate(Comparison = "FORW - SORW")) %>%
      mutate(Sector = SECTOR_LABELS[[s]], Basis = sub("elpd_", "", col))
  })
}) %>% mutate(m1 = sub(" - .*", "", Comparison),m2 = sub(".* - ", "", Comparison),
              Interpretation = case_when(z >  2 ~ paste(m1, "favored"),z < -2 ~ paste(m2, "favored"),
                                         TRUE   ~ "Indistinguishable")) %>%
  select(Sector, Comparison, Basis, Delta, SE, z, Interpretation)
cat("\n[C] 三模型两两比较 (hybrid 应与 Step 8 已有结果一致；Drift-FORW 是新增)\n")
print(as.data.frame(pair_tab %>% arrange(Comparison, Sector, Basis)))
#write_csv(pair_tab, file.path(OUT, "C_pairwise_drift_forw_sorw.csv"))

## stacking / pseudo-BMA+ 权重 (3 个 Reduced 动力学模型) + 点重抽样区间 + 加权 P_exc
set.seed(20261001)
stack_tab <- purrr::map_dfr(SECTOR_CODES, function(s) {
  w <- wide_h(s, "elpd_hybrid")
  M <- as.matrix(w[, c("Reduced_SORW", "Reduced_FORW", "Reduced_DriftFORW")])
  ws <- as.numeric(suppressMessages(loo::stacking_weights(M)))
  wb <- as.numeric(suppressMessages(loo::pseudobma_weights(M, BB = TRUE)))
  boot <- t(replicate(300, { i <- sample(nrow(M), replace = TRUE)
  as.numeric(suppressMessages(loo::stacking_weights(M[i, , drop = FALSE]))) }))
  pe <- step4_forecast$analytic_2035_p_exc %>% filter(Sector == SECTOR_LABELS[[s]])
  p_vec <- c(pe$P_exc_SORW, pe$P_exc_FORW, pe$P_exc_DRIFT)
  tibble(Sector = SECTOR_LABELS[[s]],
         w_stack_SORW = ws[1], w_stack_FORW = ws[2], w_stack_Drift = ws[3],
         w_stack_Drift_lo = quantile(boot[, 3], .1), w_stack_Drift_hi = quantile(boot[, 3], .9),
         w_BMA_SORW = wb[1], w_BMA_FORW = wb[2], w_BMA_Drift = wb[3],
         P_exc_stacked = sum(ws * p_vec),P_exc_BMA = sum(wb * p_vec), 
         P_exc_min = min(p_vec), P_exc_max = max(p_vec))
})
cat("\n[C] stacking 权重 (14 个点，权重不稳定，区间为 10-90% 点重抽样)\n")
print(as.data.frame(stack_tab))
#write_csv(stack_tab, file.path(OUT, "C_stacking_weights.csv"))

# =============================================================================
# D. LFO: LPD + CRPS + 80/95% 区间覆盖 + COVID/非COVID 拆分（对应问题 5）
#    注意：已修正 LLT 的 h 步预测方差:  Var = h*s_level^2 + sum_{j=1}^{h-1} j^2 * s_slope^2 + s_Y^2
#    (原 07 里 level 项用了 sum_{j=1}^h j^2，对 h=2 高估 (5 vs 2))。
#    origin 扩展为 2015:2021，以增加非 COVID 目标点（训练集至少 6 个点）。
#    共 5 部门 x 14 任务 x 4 模型 = 280 次拟合 (iter=2000)。
# =============================================================================
crps_sample <- function(y, x) {                       # 样本版 CRPS = E|X-y| - 0.5 E|X-X'|
  n <- length(x); xs <- sort(x); i <- seq_len(n)
  mean(abs(x - y)) - sum((2 * i - n - 1) * xs) / n^2
}
ex <- function(fit, pars) lapply(rstan::extract(fit, pars = pars), function(v) if (is.matrix(v)) v else as.numeric(v))

pp_sorw <- function(fit, h, cov_t) {
  p <- ex(fit, c("mu_trend", "s_mu", "s_Y", "beta_covid")); n <- ncol(p$mu_trend)
  list(m = p$mu_trend[, n] + h * (p$mu_trend[, n] - p$mu_trend[, n - 1]) + p$beta_covid * cov_t,
       s = sqrt(sum(seq_len(h)^2) * p$s_mu^2 + p$s_Y^2))
}
pp_forw <- function(fit, h, cov_t) {
  p <- ex(fit, c("mu_trend", "s_mu", "s_Y", "beta_covid")); n <- ncol(p$mu_trend)
  list(m = p$mu_trend[, n] + p$beta_covid * cov_t, s = sqrt(h * p$s_mu^2 + p$s_Y^2))
}
pp_drift <- function(fit, h, cov_t) {
  p <- ex(fit, c("mu_trend", "drift", "s_mu", "s_Y", "beta_covid")); n <- ncol(p$mu_trend)
  list(m = p$mu_trend[, n] + h * p$drift + p$beta_covid * cov_t, s = sqrt(h * p$s_mu^2 + p$s_Y^2))
}
pp_llt <- function(fit, h, cov_t) {
  p <- ex(fit, c("level", "slope", "s_level", "s_slope", "s_Y", "beta_covid"))
  w_slope <- if (h <= 1) 0 else sum(seq_len(h - 1)^2)
  list(m = p$level[, ncol(p$level)] + h * p$slope[, ncol(p$slope)] + p$beta_covid * cov_t,
       s = sqrt(h * p$s_level^2 + w_slope * p$s_slope^2 + p$s_Y^2))      # <- 修正处
}
score_one <- function(pp, y_obs, seed) {
  set.seed(seed)
  ys <- rnorm(length(pp$m), pp$m, pp$s)               # 每个后验抽样一个预测抽样 = 混合预测分布的样本
  q <- quantile(ys, c(.025, .10, .90, .975), names = FALSE)
  tibble(LPD = log_mean_exp(dnorm(y_obs, pp$m, pp$s, log = TRUE)), CRPS = crps_sample(y_obs, ys),
         In80 = y_obs >= q[2] & y_obs <= q[3], In95 = y_obs >= q[1] & y_obs <= q[4], Width95 = q[4] - q[1])
}
fit_llt_prefix <- function(Y_tr, seed) {
  sampling(MODELS$Reduced_LLT,
           data = list(T = length(Y_tr), Y = Y_tr, covid_dummy = as.numeric(COVID_PRIMARY[seq_along(Y_tr)]),
                       prior_alpha = CFG$prior_alpha, prior_beta = CFG$prior_beta),
           iter = 2000, seed = seed,
           control = list(adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth), refresh = 0)
}

lfo_origins2 <- 2015:2021; lfo_h <- 1:2
lfo_res <- list(); cnt <- 0
for (sector in SECTOR_CODES) {
  Y_full <- as.numeric(D[[sector]])
  for (o in lfo_origins2) {
    for (h in lfo_h) {
      ty <- o + h
      if (!(ty %in% CFG$years)) next
      cnt <- cnt + 1
      Y_tr <- Y_full[CFG$years <= o]
      cov_t <- as.numeric(ty %in% CFG$covid_years)
      y_obs <- log(Y_full[CFG$years == ty])
      pps <- list(
        SORW      = pp_sorw (fit_prefix_model(MODELS$Reduced_SORW,      Y_tr, COVID_PRIMARY, 21000 + cnt), h, cov_t),
        FORW      = pp_forw (fit_prefix_model(MODELS$Reduced_FORW,      Y_tr, COVID_PRIMARY, 22000 + cnt), h, cov_t),
        DriftFORW = pp_drift(fit_prefix_model(MODELS$Reduced_DriftFORW, Y_tr, COVID_PRIMARY, 23000 + cnt), h, cov_t),
        LLT       = pp_llt  (fit_llt_prefix(Y_tr, 24000 + cnt), h, cov_t))
      for (m in names(pps)) {
        lfo_res[[length(lfo_res) + 1]] <- bind_cols(
          tibble(Sector = SECTOR_LABELS[[sector]], Origin = o, TargetYear = ty, Horizon = h,
                 COVID_Target = cov_t == 1, Model = m),
          score_one(pps[[m]], y_obs, 99000 + cnt))
      }
      cat(sprintf("LFO %s origin %d -> %d done (%d)\n", sector, o, ty, cnt))
    }
  }
}
lfo2 <- bind_rows(lfo_res) %>% mutate(Group = case_when(!COVID_Target ~ "Non-COVID target",
                                                        Origin >= 2020 ~ "COVID target (beta identified)",
                                                        TRUE ~ "COVID target (beta prior-driven)"))
#write_csv(lfo2, file.path(OUT, "D_LFO_scores_long.csv"))
lfo_summary2 <- bind_rows(mutate(lfo2, Group = "All targets"), lfo2) %>% group_by(Group, Model) %>%
  summarise(N = n(), Mean_LPD = mean(LPD), Mean_CRPS = mean(CRPS),
            Cover80 = mean(In80), Cover95 = mean(In95), Mean_Width95 = mean(Width95), .groups = "drop") %>%
  arrange(Group, Mean_CRPS)
cat("\n[D] LFO 汇总 (CRPS 越小越好; 覆盖率名义值 0.80/0.95)\n"); print(as.data.frame(lfo_summary2))
lfo_by_sector <- bind_rows(mutate(lfo2, Group = "All targets"), filter(lfo2, Group == "Non-COVID target")) %>%
  group_by(Group, Sector, Model) %>% summarise(N = n(), Mean_LPD = mean(LPD), Mean_CRPS = mean(CRPS), .groups = "drop")
print(as.data.frame(lfo_by_sector))
#write_csv(lfo_summary2, file.path(OUT, "D_LFO_summary.csv")); write_csv(lfo_by_sector, file.path(OUT, "D_LFO_by_sector.csv"))

# =============================================================================
# E. 先验比较（对应问题 7）：旧 Gamma vs 现用 LogNormal，不需要 Stan
#    old_a/old_b 改成旧模型实际使用的 Gamma(shape, rate)
# =============================================================================
old_a <- 1; old_b <- 2
qs <- c(.025, .5, .975)
prior_tab <- bind_rows(
  tibble(Param = "sigma_mu", Prior = "old Gamma",  Q = qs, Value = qgamma(qs, old_a, old_b)),
  tibble(Param = "sigma_mu", Prior = "LogNormal",  Q = qs, Value = qlnorm(qs, -2.0, 0.7)),
  tibble(Param = "sigma_Y",  Prior = "old Gamma",  Q = qs, Value = qgamma(qs, old_a, old_b)),
  tibble(Param = "sigma_Y",  Prior = "LogNormal",  Q = qs, Value = qlnorm(qs, -2.3, 0.6))) %>%
  pivot_wider(names_from = Q, values_from = Value, names_prefix = "q")
mass_tab <- tibble(
  Param = c("sigma_mu", "sigma_mu", "sigma_Y", "sigma_Y"),
  Prior = c("old Gamma", "LogNormal", "old Gamma", "LogNormal"),
  P_lt_0.01 = c(pgamma(.01, old_a, old_b), plnorm(.01, -2.0, .7), pgamma(.01, old_a, old_b), plnorm(.01, -2.3, .6)),
  P_lt_0.02 = c(pgamma(.02, old_a, old_b), plnorm(.02, -2.0, .7), pgamma(.02, old_a, old_b), plnorm(.02, -2.3, .6)))
cat("\n[E] 先验分位数\n"); print(as.data.frame(prior_tab)); print(as.data.frame(mass_tab))
set.seed(1); N <- 2e5
R_old <- local({ a <- rgamma(N, old_a, old_b); b <- rgamma(N, old_a, old_b); a^2 / (a^2 + b^2) })
cat("old Gamma 下 R_pers prior 2.5/50/97.5% =", round(quantile(R_old, qs), 3), "\n")

post_sig <- purrr::map_dfr(SECTOR_CODES, function(s) {
  p <- rstan::extract(step2_results$fits[[s]]$Reduced_SORW, pars = c("s_mu", "s_Y"))
  tibble(Sector = SECTOR_LABELS[[s]], s_mu = as.numeric(p$s_mu), s_Y = as.numeric(p$s_Y))
}) %>% pivot_longer(c(s_mu, s_Y), names_to = "Par", values_to = "sigma")
grid <- expand_grid(sigma = exp(seq(log(0.003), log(1.5), length.out = 400)), Par = c("s_mu", "s_Y")) %>%
  mutate(LogNormal = sigma * dlnorm(sigma, if_else(Par == "s_mu", -2.0, -2.3), if_else(Par == "s_mu", 0.7, 0.6)),
         Gamma_old = sigma * dgamma(sigma, old_a, old_b)) %>%      # 乘 sigma = 变换到 log(sigma) 尺度的密度
  pivot_longer(c(LogNormal, Gamma_old), names_to = "Prior", values_to = "dens")
p_prior <- ggplot() + geom_density(data = post_sig, aes(log(sigma), colour = Sector), linewidth = 0.5) +
  geom_line(data = grid, aes(log(sigma), dens, linetype = Prior), linewidth = 0.9) +
  facet_wrap(~Par) + labs(x = "log(sigma)", y = "density", title = "Prior vs posterior (Reduced SORW)") + theme_bw()
ggsave(file.path(OUT, "E_prior_vs_posterior_sigma.png"), p_prior, width = 9, height = 4, dpi = 300)

# =============================================================================
# F. 请把下面这些输出也发给我（问题 6 / 干预图）
# =============================================================================
cat("\n[F] 干预基线与 shock-only 表\n")
print(filter(step6_table, Scenario == "Baseline"), n = Inf)
print(shock_only_table, n = Inf)
cat("\n[F] 协变量共线性审计\n")
if (exists("collinearity_audit")) print(as.data.frame(collinearity_audit))
if (exists("cor_X_detrended"))   print(round(cor_X_detrended, 3))
cat("\nsupplementary_checks.R complete\n")


# supplementary_checks_3.R
# Run in the session that holds the objects of summaryAllStepCode.R.
# Needs: step2_results, simulation_df, SECTOR_CODES, SECTOR_LABELS, D, COVID_PRIMARY, OUTPUT_DIR, fit_sorw_grid
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })
OUT <- file.path(OUTPUT_DIR, "supp"); dir.create(OUT, showWarnings = FALSE)

# --- A. Log width W_h of the 95% predictive interval (SORW, FORW, DriftFORW) ---------------
set.seed(2026)
width_one <- function(fit, model, h_max = 12, reps = 3) {
  p <- rstan::extract(fit, pars = c("mu_trend", "s_mu", "s_Y", if (model == "DriftFORW") "drift"))
  n <- ncol(p$mu_trend); muT <- p$mu_trend[, n]; muT1 <- p$mu_trend[, n - 1]
  purrr::map_dfr(seq_len(h_max), function(h) {
    m <- switch(model, SORW = muT + h * (muT - muT1), FORW = muT, DriftFORW = muT + h * p$drift)
    v <- switch(model, SORW = sum(seq_len(h)^2) * p$s_mu^2, FORW = h * p$s_mu^2, DriftFORW = h * p$s_mu^2) + p$s_Y^2
    y <- rnorm(length(m) * reps, rep(m, reps), sqrt(rep(v, reps)))
    q <- quantile(y, c(.025, .975), names = FALSE)
    tibble(Year = 2023 + h, Model = model, W = q[2] - q[1])
  })
}
W_tab <- purrr::map_dfr(SECTOR_CODES, function(s) {
  f <- step2_results$fits[[s]]
  bind_rows(width_one(f$Reduced_SORW, "SORW"), width_one(f$Reduced_FORW, "FORW"),
            width_one(f$Reduced_DriftFORW, "DriftFORW")) %>% mutate(Sector = SECTOR_LABELS[[s]])
})
#write_csv(W_tab, file.path(OUT, "S_logwidth_W.csv"))
print(W_tab %>% filter(Year %in% c(2024, 2027, 2030, 2035)) %>%
        pivot_wider(names_from = Model, values_from = W) %>% arrange(Sector, Year), n = Inf)
p_w <- ggplot(W_tab, aes(Year, W, colour = Model)) + geom_line(linewidth = 0.8) +
  facet_wrap(~Sector, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Log width of the 95% predictive interval") + theme_bw(base_size = 12)
ggsave(file.path(OUT, "Figure_S_LongHorizon_Uncertainty.png"), p_w, width = 10, height = 6, dpi = 300)

# --- B. Simulation: valid fits, agreement of fitted and true probability ------------------
sim_sum <- simulation_df %>% group_by(Sector, DGP) %>%
  summarise(N_total = n(), N_valid = sum(Fit_OK),
            True = mean(True_P_exc[Fit_OK]), Fitted = mean(Fitted_P_exc[Fit_OK]),
            Error = mean(P_exc_Error[Fit_OK]), RMSE = sqrt(mean(P_exc_Error[Fit_OK]^2)),
            Cor = cor(True_P_exc[Fit_OK], Fitted_P_exc[Fit_OK]),
            Share_within_0.1 = mean(abs(P_exc_Error[Fit_OK]) < 0.1), .groups = "drop")
print(as.data.frame(sim_sum))
#write_csv(sim_sum, file.path(OUT, "S_simulation_summary.csv"))

# error by bin of true probability, SORW process only (shows pull toward the middle at the extremes)
bins <- simulation_df %>% filter(Fit_OK, DGP == "SORW_correct") %>%
  mutate(Bin = cut(True_P_exc, c(0, .05, .25, .75, .95, 1), include.lowest = TRUE)) %>%
  group_by(Bin) %>% summarise(N = n(), Mean_true = mean(True_P_exc), Mean_fitted = mean(Fitted_P_exc))
print(as.data.frame(bins))

p_sim <- simulation_df %>% filter(Fit_OK) %>%
  ggplot(aes(True_P_exc, Fitted_P_exc)) + geom_abline(linetype = 2) + geom_point(alpha = 0.3, size = 0.8) +
  facet_grid(Sector ~ DGP) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "True probability", y = "Fitted probability") + theme_bw(base_size = 11)
ggsave(file.path(OUT, "Figure_S_Simulation_TrueVsFitted.png"), p_sim, width = 11, height = 5, dpi = 300)

# --- C. Grid solution against the Stan result for the empirical series --------------------
grid_check <- purrr::map_dfr(c("coal", "oil"), function(s) {
  tibble(Sector = SECTOR_LABELS[[s]],
         Grid_Pexc = fit_sorw_grid(as.numeric(D[[s]]), covid = as.numeric(COVID_PRIMARY))$P_exc)
}) %>% left_join(step4_forecast$analytic_2035_p_exc %>% select(Sector, Stan_Pexc = P_exc_SORW), by = "Sector")
print(as.data.frame(grid_check))

# =============================================================================
# 重画 Figure S: 阻尼对超标概率的边际效应 (5 载体，相对 Shock-only)
# =============================================================================
stopifnot(exists("step6_table"))

# 1. 提取 shock-only (d = 0) 基准并计算相对差值
shock_base <- step6_table %>% filter(Damping == 0) %>% select(Sector, Reduction_pct, P_exc_base = P_exc)
fig_damping_data <- step6_table %>% filter(Damping > 0) %>% # 仅保留 d = 0.25, 0.50, 0.75
  left_join(shock_base, by = c("Sector", "Reduction_pct")) %>%
  mutate(Delta_Pexc = P_exc - P_exc_base, Sector = factor(Sector, levels = SECTOR_ORDER_LABELS),
         Damping_Label = factor(Damping,levels = c(0.25, 0.50, 0.75),
                                labels = c("d = 0.25", "d = 0.50", "d = 0.75")))

# 2. 绘制分面折线图
p_damping_si <- ggplot(fig_damping_data, aes(x = Reduction_pct, y = Delta_Pexc, color = Damping_Label, group = Damping_Label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_line(linewidth = 0.85) + geom_point(size = 1.6) + facet_wrap(~ Sector, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = c("d = 0.25" = "#4F81BD", "d = 0.50" = "#2CA02C", "d = 0.75" = "#C00000")) +
  labs(x = "Permanent Consumption Reduction Level from 2025 (%)",
       y = expression(Delta * P[exc] ~ "(Relative to Shock-Only Baseline)"),color = "Innovation Damping (d)") +
  theme_bw(base_size = 11) + theme(legend.position = "top", panel.grid.minor = element_blank(),
                                   strip.background = element_rect(fill = "gray95"),
                                   strip.text = element_text(face = "bold"))
# 3. 输出保存
ggsave(file.path(OUT, "Figure_S_Intervention_Damping.png"), p_damping_si, width = 9.5, height = 5.5, dpi = 300)

cat(">>> Figure S (Intervention Damping) generated successfully.\n")

# supplementary_checks_4.R
# Run in the session that holds the objects of summaryAllStepCode.R.
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })
OUT <- file.path(OUTPUT_DIR, "supp"); dir.create(OUT, showWarnings = FALSE)
q <- function(x) sprintf("%.3f (%.3f to %.3f)", median(x), quantile(x, .025), quantile(x, .975))

# --- A. Posterior summaries of the model parameters (Reduced models) ----------------------
par_tab <- purrr::map_dfr(SECTOR_CODES, function(s) {
  purrr::map_dfr(c("Reduced_SORW", "Reduced_FORW", "Reduced_DriftFORW"), function(m) {
    f <- step2_results$fits[[s]][[m]]
    p <- rstan::extract(f, pars = c("s_mu", "s_Y", "beta_covid", "mu_trend", if (m == "Reduced_DriftFORW") "drift"))
    n <- ncol(p$mu_trend)
    tibble(Carrier = SECTOR_LABELS[[s]], Model = sub("Reduced_", "", m),
           sigma_mu = q(p$s_mu), sigma_Y = q(p$s_Y), beta_COVID = q(p$beta_covid),
           Drift = if (m == "Reduced_DriftFORW") q(p$drift) else NA_character_,
           Slope_at_T = if (m == "Reduced_SORW") q(p$mu_trend[, n] - p$mu_trend[, n - 1]) else NA_character_)
  })
})
print(as.data.frame(par_tab))
#write_csv(par_tab, file.path(OUT, "S_parameter_summaries.csv"))

# --- B. Exceedance probability under alternative anchors of the target --------------------
# (1) observed 2023 value (current), (2) latent level at 2023 (draw-specific), (3) mean of the log observations 2021-2023
pexc_anchor <- purrr::map_dfr(SECTOR_CODES, function(s) {
  y <- log(as.numeric(D[[s]])); T <- length(y); lr <- log(CFG$target_ratio); h <- 12
  purrr::map_dfr(c("SORW", "FORW", "DriftFORW"), function(m) {
    f <- step2_results$fits[[s]][[paste0("Reduced_", m)]]
    p <- rstan::extract(f, pars = c("mu_trend", "s_mu", "s_Y", if (m == "DriftFORW") "drift"))
    muT <- p$mu_trend[, T]; muT1 <- p$mu_trend[, T - 1]
    shift <- switch(m, SORW = h * (muT - muT1), FORW = 0 * muT, DriftFORW = h * p$drift)   # m_h - mu_T
    v <- switch(m, SORW = sum(seq_len(h)^2) * p$s_mu^2, FORW = h * p$s_mu^2, DriftFORW = h * p$s_mu^2) + p$s_Y^2
    mh <- muT + shift
    P <- function(anchor) mean(1 - pnorm((lr + anchor - mh) / sqrt(v)))
    tibble(Carrier = SECTOR_LABELS[[s]], Model = m,
           Observed_2023 = P(y[T]), Latent_2023 = P(muT), Mean_2021_2023 = P(mean(y[(T - 2):T])))
  })
})
print(as.data.frame(pexc_anchor))
#write_csv(pexc_anchor, file.path(OUT, "S_pexc_anchor_sensitivity.csv"))

# --- C. Exceedance probability against the target ratio (SORW, FORW) at h = 12 -------------
ratio_tab <- step4_forecast$target_ratio_results %>%
  filter(Horizon == 12, round(Target_Ratio, 2) %in% c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2)) %>%
  mutate(Target_Ratio = round(Target_Ratio, 2)) %>%
  pivot_wider(names_from = Target_Ratio, values_from = P_exc) %>% arrange(Sector, Model)
print(as.data.frame(ratio_tab))
#write_csv(ratio_tab, file.path(OUT, "S_pexc_target_ratio.csv"))