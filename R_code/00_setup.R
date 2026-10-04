# =============================================================================
# 00_setup.R
# Global initialization: libraries, paths, CFG, SECTORS, data, covariates,
# shared helpers, and Stan models.
# (全局初始化：库、路径、CFG、SECTORS、数据、协变量、共享 helper、Stan 模型)
# Every downstream module sources this file once; nothing else re-declares
# (上游所有模块只 source 一次这个文件，不再重复 library()/read_csv()/stan_model())
# =============================================================================
suppressPackageStartupMessages({
  library(rstan)
  library(loo)
  library(tidyverse)
  library(matrixStats)
  library(forcats)
  library(purrr)
  library(tibble)
  library(dplyr)
})
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# --- Paths (路径) -------------------------------------------------------------------
STAN_DIR   <- "~/Jo/CSUC/PaperSubmition/ZRX/EMS_ Manuscript/Analysis/"
OUTPUT_DIR <- file.path(STAN_DIR, "model_outputs")
FIG_DIR    <- file.path(STAN_DIR, "figs")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,    recursive = TRUE, showWarnings = FALSE)

# --- Global configuration (全局配置) ---------------------------------------------------------------
CFG <- list(years= 2010:2023,future_years= 2024:2035,target_ratio= 0.80,intervention_year = 2025,
            covid_years= c(2020, 2021, 2022),seed= 1234,iter= 4000,warmup= 2000,chains= 4,
            adapt_delta= 0.99,max_treedepth= 14,prior_alpha= 1,prior_beta= 2,lasso_alpha= 30,lasso_beta= 30)

# --- Sector codes and labels (Sector 代码与标签)------------------------------------------
SECTOR_CODES  <- c("coal", "oil", "NGas", "thermal", "electric")
SECTOR_LABELS <- c(coal= "Coal", oil= "Oil Products", NGas= "Natural Gas",
                   thermal= "Thermal (Heat)", electric= "Electricity")
SECTOR_ORDER_LABELS <- c("Coal", "Oil Products", "Natural Gas", "Thermal (Heat)", "Electricity")
# --- Data loading (数据加载) ---------------------------------------------------------------
DIPAT <- read_csv("~/Jo/Hirosaki Univ./zhaoruixi_col/data2025/EP/appendix/tableA4.csv",show_col_types = FALSE)
D <- read.csv("~/Jo/Hirosaki Univ./zhaoruixi_col/data2025/EP/appendix/tableA5.csv")
stopifnot(nrow(DIPAT) == length(CFG$years))
stopifnot(nrow(D)     == length(CFG$years))
stopifnot(all(SECTOR_CODES %in% colnames(D)))
# --- Covariates (协变量) -----------------------------------------------------------------
X_LOG <- DIPAT %>% select(P, A, Tstr, Ttech) %>% mutate(across(everything(), log)) %>% as.matrix()
stopifnot(nrow(X_LOG) == length(CFG$years), ncol(X_LOG) == 4, all(is.finite(X_LOG)))
time_index <- seq_len(nrow(X_LOG))
X_FULL <- as.matrix(as.data.frame(lapply(as.data.frame(X_LOG), function(x) resid(lm(x ~ time_index)))))
colnames(X_FULL) <- colnames(X_LOG)
# --- Primary COVID coding (COVID 主编码) -----------------------------------------------------------
COVID_PRIMARY <- ifelse(CFG$years %in% CFG$covid_years, 1, 0)
# =============================================================================
# Shared helpers (共享 helper)
# =============================================================================
# --- Stan data constructors (Stan data 构造函数) -----------------------------------------------------
make_full_data <- function(Y, X = X_FULL, covid_dummy = COVID_PRIMARY, prior = CFG) {
  stopifnot(length(Y) == length(CFG$years), all(is.finite(Y)), all(Y > 0))
  list(T = length(Y), K = ncol(X), X = X, Y = as.numeric(Y), covid_dummy = as.numeric(covid_dummy),
       prior_alpha = prior$prior_alpha, prior_beta = prior$prior_beta,
       lasso_alpha = prior$lasso_alpha, lasso_beta = prior$lasso_beta)
}
make_reduced_data <- function(Y, covid_dummy = COVID_PRIMARY, prior = CFG) {
  stopifnot(length(Y) == length(CFG$years), all(is.finite(Y)), all(Y > 0))
  list(T = length(Y), Y = as.numeric(Y),covid_dummy = as.numeric(covid_dummy),
       prior_alpha = prior$prior_alpha, prior_beta = prior$prior_beta)
}
# --- Diagnostic function (诊断函数)-------------------------------------------------
diagnose_stan_fit <- function(fit, model_name = NA_character_, sector = NA_character_) {
  summ <- summary(fit)$summary
  rhat <- summ[, "Rhat"]
  neff <- if ("n_eff" %in% colnames(summ)) summ[, "n_eff"] else rep(NA_real_, length(rhat))
  
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences    <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= CFG$max_treedepth), numeric(1)))
  energy <- lapply(sp, function(x) x[, "energy__"])
  bfmi   <- vapply(energy, function(e) mean(diff(e)^2) / var(e), numeric(1))
  total_draws <- sum(vapply(sp, nrow, numeric(1)))
  div_rate    <- divergences / total_draws
  
  rhat_ok <- all(rhat < 1.01, na.rm = TRUE)
  bfmi_ok <- all(bfmi > 0.2, na.rm = TRUE)
  div_ok  <- div_rate < 0.01
  
  tibble(Sector= sector, Model= model_name, Max_Rhat= max(rhat, na.rm = TRUE), Min_n_eff= min(neff, na.rm = TRUE),
         Divergences= divergences, Max_Treedepth_Hits = treedepth_hits, Min_BFMI= min(bfmi, na.rm = TRUE),
         Rhat_OK= rhat_ok, Divergence_OK= div_ok,
         Treedepth_OK= treedepth_hits == 0, BFMI_OK= bfmi_ok)
}

RETRY_LOG <- list()
fit_retry <- function(model, data, seed, label = "", iter = CFG$iter, delta_grid = c(CFG$adapt_delta, 0.995, 0.999),
                      max_treedepth = CFG$max_treedepth) {
  best <- NULL; best_score <- Inf; best_i <- NA_integer_; passed <- FALSE
  for (i in seq_along(delta_grid)) {
    fit <- tryCatch(
      sampling(model, data = data, iter = iter, seed = seed + 500000L * (i - 1L),
               control = list(adapt_delta = delta_grid[i], max_treedepth = max_treedepth), refresh = 0),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    d <- diagnose_stan_fit(fit)
    ok <- d$Rhat_OK && d$Divergence_OK && d$BFMI_OK
    if (ok) { best <- fit; best_i <- i; passed <- TRUE; break }
    score <- (!d$Rhat_OK) * 1e6 + (!d$BFMI_OK) * 1e4 + d$Divergences
    if (score < best_score) { best <- fit; best_score <- score; best_i <- i }
  }
  RETRY_LOG[[length(RETRY_LOG) + 1]] <<- tibble(Label = label, Attempt_Used = best_i, Passed = passed)
  if (!passed) warning("No attempt passed diagnostics: ", label)
  best
}
# --- Pointwise ELPD comparison (比较) ----------------------------------------------------
compare_pointwise <- function(loo_1, loo_2, model_1, model_2) {
  pw1 <- as.data.frame(loo_1$pointwise)
  pw2 <- as.data.frame(loo_2$pointwise)
  delta_i  <- pw1$elpd_loo - pw2$elpd_loo
  delta    <- sum(delta_i)
  se_delta <- sqrt(length(delta_i) * var(delta_i))
  tibble(Model_1 = model_1, Model_2 = model_2,
         Delta_ELPD = delta, SE_Delta = se_delta, Delta_over_SE = delta / se_delta)
}

# --- Evidence classification threshold (分类阈值) ---------------------------------------------------------------
classify_evidence <- function(z, pos_label, neg_label, threshold = 2) {
  dplyr::case_when(is.na(z) ~ "Not available",z >  threshold ~ pos_label,
                   z < -threshold ~ neg_label,TRUE ~ "No clear predictive preference")
}

# --- Posterior summary (后验总结) ---------------------------------------------------------------
summarize_posterior <- function(x, prob = 0.95) {
  qs <- quantile(x, probs = c((1 - prob)/2, 0.5, 1 - (1 - prob)/2), na.rm = TRUE)
  data.frame(Mean = mean(x, na.rm = TRUE), Median = qs[2], Lower = qs[1], Upper = qs[3])
}

# ---  log-mean-exp ------------------------------------------------------
log_mean_exp <- function(x) {
  if (length(x) == 0) stop("Empty input to log_mean_exp")
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# =============================================================================
# Stan model compilation (Stan 模型编译)
# v2 = primary analysis models; v3 = exact-LOO models (with use_obs mask)
# v2 = 主分析模型；v3 = exact-LOO 专用（含 use_obs 掩码）
# Both versions compiled here to avoid recompilation in 08_exact_LOO.R
# (两个版本并行编译，避免 08_exact_LOO.R 重新编译)
# =============================================================================
MODELS <- list(
  Full_SORW            = stan_model(file.path(STAN_DIR, "SSM_Full_SORW_v2.stan")),
  Reduced_SORW         = stan_model(file.path(STAN_DIR, "SSM_Reduced_SORW_v2.stan")),
  Reduced_FORW         = stan_model(file.path(STAN_DIR, "SSM_Reduced_FORW_v2.stan")),
  Reduced_DriftFORW    = stan_model(file.path(STAN_DIR, "SSM_Reduced_DriftFORW_v1.stan")),
  Reduced_LLT          = stan_model(file.path(STAN_DIR, "SSM_Reduced_LLT_v1.stan")),
  Full_SORW_v3         = stan_model(file.path(STAN_DIR, "SSM_Full_SORW_v3.stan")),
  Reduced_SORW_v3      = stan_model(file.path(STAN_DIR, "SSM_Reduced_SORW_v3.stan")),
  Reduced_FORW_v3      = stan_model(file.path(STAN_DIR, "SSM_Reduced_FORW_v3.stan")),
  Reduced_DriftFORW_v3 = stan_model(file.path(STAN_DIR, "SSM_Reduced_DriftFORW_v3.stan"))
)

cat("=== 00_setup.R complete ===\n")