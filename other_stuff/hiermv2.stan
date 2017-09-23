data {
  int<lower=0> N;
  int<lower=1> K;
  int<lower=1> J_station;
  int<lower=1> J_wkday;
  int<lower=1, upper=J_station> jj_station[N];
  int<lower=1, upper=J_wkday> jj_wkday[N];
  matrix[N, K] x_wkday;
  vector[N] log_count;
  real<lower=0> scale_wip;
  real<lower=0> scale_sig;
  int<lower = 0, upper = 1> run_estimation;
}
parameters {
  // z refers to standardized version of the variable
  real a_z;
  real<lower=0> sig_a_station;
  real mu_a_station;
  vector[J_station] a_station_z;
  vector[K] mu_wkday;
  vector<lower=0>[K] tau_wkday;
  cholesky_factor_corr[K] L_Omega_wkday;
  vector[K] b_wkday_z[J_wkday];
  real<lower=0> sig_log_count;
}
transformed parameters {
  vector[J_station] a_station;
  matrix[K, K] L;
  vector[K] b_wkday[J_wkday];
  real a;

  a_station = mu_a_station + sig_a_station * a_station_z;

  L = diag_pre_multiply(tau_wkday, L_Omega_wkday);
  for (j in 1:J_wkday)
    b_wkday[j] = mu_wkday + L * b_wkday_z[j];

  a = mean(log_count) + 2*sd(log_count) * a_z;
}
model {
  mu_a_station ~ normal(0, 0.01);
  sig_a_station ~ cauchy(0, scale_sig);
  a_station_z ~ normal(0, 1);

  mu_wkday[1] ~ normal(0, 0.01);
  for (k in 2:K)
    mu_wkday[k] ~ normal(0, scale_wip);
  tau_wkday ~ cauchy(0, scale_sig);
  L_Omega_wkday ~ lkj_corr_cholesky(2);
  for (j in 1:J_wkday)
    b_wkday_z[j] ~ normal(0, 1);

  a_z ~ normal(0, 1);
  sig_log_count ~ cauchy(0, scale_sig);
  if (run_estimation == 1) {
    vector[N] wk_i;
    for (n in 1:N)
      wk_i[n] = x_wkday[n] * b_wkday[jj_wkday[n]];

    log_count ~ normal(a + a_station[jj_station] + wk_i, sig_log_count);
  }
}
generated quantities {
  vector[N] y_sim;

  for (n in 1:N)
    y_sim[n] = normal_rng(a + a_station[jj_station[n]] +
                            x_wkday[n] * b_wkday[jj_wkday[n]], sig_log_count);
}
