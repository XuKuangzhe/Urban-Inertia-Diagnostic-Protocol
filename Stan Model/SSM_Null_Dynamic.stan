data {
  int<lower=1> T;             // Time series length
  vector<lower=0>[T] Y;       // Dependent variable
  vector[T] covid_dummy;
  real<lower=0> prior_alpha;  // Hyperparameter for Gamma prior
  real<lower=0> prior_beta;   // Hyperparameter for Gamma prior
}

transformed data {
  vector[T] log_Y = log(Y); 
}

parameters {
  vector[T] mu_trend;         // Latent inertial trend
  real beta_covid;
  real<lower=0> s_mu;         // Trend noise scale
  real<lower=0> s_Y;          // Observation noise scale
}

transformed parameters {
  vector[T] mu_total = mu_trend + beta_covid * covid_dummy; // No exogenous regressors
}

model {
  // 1. Variance Priors
  s_mu ~ gamma(prior_alpha, prior_beta);
  s_Y  ~ gamma(prior_alpha, prior_beta);
  
  // We expect a negative shock, but let data decide. 
  // Normal(0, 1) allows for reasonable magnitude on log scale.
  beta_covid ~ normal(0, 1);
  
  // 2. Second-Order Random Walk (Inertial Trend)
  mu_trend[1] ~ normal(log_Y[1], 1);
  mu_trend[2] ~ normal(mu_trend[1], s_mu);
  for (t in 3:T) {
    mu_trend[t] ~ normal(2*mu_trend[t-1] - mu_trend[t-2], s_mu);
  }
  
  // 3. Observation Likelihood
  log_Y ~ normal(mu_total, s_Y);
}

generated quantities {
  vector[T] log_lik;
  for(n in 1:T)
   log_lik[n] = normal_lpdf(log_Y[n] | mu_total[n], s_Y);
}
