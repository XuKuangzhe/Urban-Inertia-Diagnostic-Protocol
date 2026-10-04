// SSM_Full_SORW_v3.stan
data {
  int<lower=3> T;                  // Time-series length
  int<lower=1> K;                  // Number of socioeconomic predictors
  matrix[T, K] X;                  // Log-transformed socioeconomic predictors
  vector<lower=0>[T] Y;            // Energy consumption
  vector[T] covid_dummy;           // COVID intervention indicator

  // 1 = observation contributes to likelihood
  // 0 = observation is held out for exact LOO
  int<lower=0, upper=1> use_obs[T];

  real<lower=0> prior_alpha;       // Gamma prior shape
  real<lower=0> prior_beta;        // Gamma prior rate

  // Hyperprior for global LASSO scale
  real<lower=0> lasso_alpha;
  real<lower=0> lasso_beta;
}

transformed data {
  vector[T] log_Y = log(Y);
  real prior_scale_ratio = (prior_alpha / prior_beta) / 0.5;
  real log_prior_shift   = log(prior_scale_ratio);
}

parameters {
  //vector[K] beta; // Socioeconomic regression coefficients
  vector[K] beta_raw;
  real beta_covid; // COVID intervention effect
  // Latent SORW initial conditions
  real mu1;
  real delta1;
  vector[T - 2] z_mu; // Non-centered SORW innovations
  // Scale parameters
  real<lower=0> s_mu;
  real<lower=0> s_Y;
  real<lower=0> S_beta; // Global Bayesian LASSO scale
}

transformed parameters {
  vector[K] beta = S_beta * beta_raw;
  vector[T] mu_trend;
  vector[T] mu_total;
  vector[T - 1] slope;

  slope[1] = delta1;
  for (t in 2:(T - 1))
    slope[t] = slope[t - 1] + s_mu * z_mu[t - 1];

  mu_trend[1] = mu1;
  for (t in 2:T)
    mu_trend[t] = mu_trend[t - 1] + slope[t - 1];

  mu_total = mu_trend + covid_dummy * beta_covid + X * beta;
}

model {
  // Variance priors
  s_Y  ~ lognormal(-2.3 + log_prior_shift, 0.6);
  s_mu ~ lognormal(-2.0 + log_prior_shift, 0.7);
  // Hierarchical Bayesian LASSO
  S_beta ~ gamma(lasso_alpha, lasso_beta);
  beta_raw ~ double_exponential(0, 1);
  beta_covid ~ normal(0, 1); // COVID intervention
  mu1 ~ normal(0, 5); // Initial latent state
  delta1 ~ normal(0, 1); // Initial slope / first difference
  z_mu ~ normal(0, 1); // Non-centered SORW innovations
  // Observation model
  // Held-out observations do NOT contribute to the posterior.
  for (t in 1:T) {
    if (use_obs[t] == 1) {
      log_Y[t] ~ normal(mu_total[t], s_Y);
    }
  }
}

generated quantities {
  vector[T] log_lik;

  for (t in 1:T)
    log_lik[t] = normal_lpdf(log_Y[t] | mu_total[t], s_Y);
}
