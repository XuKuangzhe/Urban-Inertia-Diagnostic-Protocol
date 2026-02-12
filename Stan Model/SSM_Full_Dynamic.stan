data {
  int<lower=1> T;             // Time series length
  int<lower=1> K;             // Number of predictors
  matrix[T, K] X;             // Predictor matrix (log-transformed)
  vector<lower=0>[T] Y;       // Dependent variable (Energy consumption)
  real<lower=0> prior_alpha;  // Hyperparameter for Gamma prior (shape)
  real<lower=0> prior_beta;   // Hyperparameter for Gamma prior (rate)
}

transformed data {
  vector[T] log_Y = log(Y);   // Log-transform outcome
}

parameters {
  vector[K] beta;             // Regression coefficients (Elasticities)
  vector[T] mu_trend;         // Latent inertial trend
  real<lower=0> s_mu;         // Trend noise scale
  real<lower=0> s_Y;          // Observation noise scale
  real<lower=0> S_beta;       // Shrinkage parameter for Bayesian LASSO
}

transformed parameters {
  vector[T] mu_total = mu_trend + X * beta;  // Total estimated level
}

model {
  // 1. Variance Priors (Regularization)
  s_mu ~ gamma(prior_alpha, prior_beta);
  s_Y  ~ gamma(prior_alpha, prior_beta);
  
  // 2. Bayesian LASSO Hierarchical Prior
  S_beta ~ gamma(30, 30);                // Strong anchor around 1
  beta ~ double_exponential(0, S_beta);  // Laplace prior for shrinkage
  
  // 3. Second-Order Random Walk (Inertial Trend)
  mu_trend[1] ~ normal(log_Y[1], 1);
  mu_trend[2] ~ normal(mu_trend[1], s_mu);
  for (t in 3:T) {
    mu_trend[t] ~ normal(2*mu_trend[t-1] - mu_trend[t-2], s_mu);
  }
  
  // 4. Observation Likelihood
  log_Y ~ normal(mu_total, s_Y);
}

generated quantities {
  vector[T] log_lik;
  for(n in 1:T)
   log_lik[n] = normal_lpdf(log_Y[n] | mu_total[n], s_Y);
}


