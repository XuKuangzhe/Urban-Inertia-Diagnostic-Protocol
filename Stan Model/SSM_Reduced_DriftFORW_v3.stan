// SSM_Reduced_DriftFORW_v3.stan  (exact-LOO version with use_obs mask)
data {
  int<lower=2> T;
  vector<lower=0>[T] Y;
  vector[T] covid_dummy;
  array[T] int<lower=0, upper=1> use_obs;
  real<lower=0> prior_alpha;
  real<lower=0> prior_beta;
}
transformed data {
  vector[T] log_Y = log(Y);
  real log_prior_shift = log((prior_alpha / prior_beta) / 0.5);
}
parameters {
  real beta_covid;
  real mu1;
  real drift;
  vector[T - 1] z_mu;
  real<lower=0> s_mu;
  real<lower=0> s_Y;
}
transformed parameters {
  vector[T] mu_trend;
  vector[T] mu_total;
  mu_trend[1] = mu1;
  for (t in 2:T)
    mu_trend[t] = mu_trend[t - 1] + drift + s_mu * z_mu[t - 1];
  mu_total = mu_trend + covid_dummy * beta_covid;
}
model {
  s_mu ~ lognormal(-2.0 + log_prior_shift, 0.7);
  s_Y  ~ lognormal(-2.3 + log_prior_shift, 0.6);
  beta_covid ~ normal(0, 1);
  mu1 ~ normal(0, 5);
  drift ~ normal(0, 1);
  z_mu ~ normal(0, 1);
  for (t in 1:T)
    if (use_obs[t] == 1)
      target += normal_lpdf(log_Y[t] | mu_total[t], s_Y);
}
generated quantities {
  vector[T] log_lik;
  for (t in 1:T)
    log_lik[t] = normal_lpdf(log_Y[t] | mu_total[t], s_Y);
}
