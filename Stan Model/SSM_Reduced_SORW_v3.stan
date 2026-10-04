// SSM_Reduced_SORW_v3.stan
data {
  int<lower=3> T;
  vector<lower=0>[T] Y;
  vector[T] covid_dummy;
  real<lower=0> prior_alpha;
  real<lower=0> prior_beta;
  int<lower=0, upper=1> use_obs[T];
}

transformed data {
  vector[T] log_Y = log(Y);
  real prior_scale_ratio = (prior_alpha / prior_beta) / 0.5;
  real log_prior_shift   = log(prior_scale_ratio);
}

parameters {
  real beta_covid; // COVID intervention effect
  // Latent SORW initial conditions
  real mu1;
  real delta1;
  vector[T - 2] z_mu; // Non-centered SORW innovations
  // Scale parameters
  real<lower=0> s_mu;
  real<lower=0> s_Y;
}

transformed parameters {
  vector[T] mu_trend;
  vector[T] mu_total;
  vector[T - 1] slope;

  slope[1] = delta1;
  for (t in 2:(T - 1))
    slope[t] = slope[t - 1] + s_mu * z_mu[t - 1];

  mu_trend[1] = mu1;
  for (t in 2:T)
    mu_trend[t] = mu_trend[t - 1] + slope[t - 1];

  mu_total = mu_trend + covid_dummy * beta_covid;
}

model {
  // Variance priors
  s_Y  ~ lognormal(-2.3 + log_prior_shift, 0.6);
  s_mu ~ lognormal(-2.0 + log_prior_shift, 0.7);
  beta_covid ~ normal(0, 1); // COVID intervention
  // Initial latent state
  //mu1 ~ normal(log_Y[1], 1);
  mu1 ~ normal(0, 5);
  delta1 ~ normal(0, 1);
  z_mu ~ normal(0, 1); // Non-centered SORW innovations
  // Observation modellog_Y ~ normal(mu_total, s_Y);
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
