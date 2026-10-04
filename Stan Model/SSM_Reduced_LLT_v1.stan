//SSM_Reduced_LLT_v1.stan
data {
  int<lower=2> T;
  vector<lower=0>[T] Y;
  vector[T] covid_dummy;
  real<lower=0> prior_alpha;
  real<lower=0> prior_beta;
}

transformed data {
  vector[T] log_Y = log(Y);
}

parameters {
  real beta_covid;
  real mu1;
  real delta1;
  vector[T - 1] z_level;
  vector[T - 1] z_slope;
  real<lower=0> s_level;
  real<lower=0> s_slope;
  real<lower=0> s_Y;
}

transformed parameters {
  vector[T] level;
  vector[T] slope;
  vector[T] mu_total;
  level[1] = mu1;
  slope[1] = delta1;

  for (t in 2:T) {
    level[t] = level[t - 1]+ slope[t - 1] + s_level * z_level[t - 1];
    slope[t] =slope[t - 1]+ s_slope * z_slope[t - 1];
  }

  mu_total =level + covid_dummy * beta_covid;
}

model {
  s_Y     ~ lognormal(-2.3, 0.6);
  s_level ~ lognormal(-2.0, 0.7);
  s_slope ~ lognormal(-2.0, 0.7);
  beta_covid ~ normal(0, 1);
  mu1 ~ normal(0, 5);
  delta1 ~ normal(0, 1);
  z_level ~ normal(0, 1);
  z_slope ~ normal(0, 1);

  log_Y ~ normal(mu_total,s_Y);
}

generated quantities {
  vector[T] log_lik;
  for (t in 1:T)
    log_lik[t] =normal_lpdf(log_Y[t] |mu_total[t],s_Y);
}
