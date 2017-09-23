data {
  int<lower=0> N;
  int<lower=1> J_station;
  int<lower=1> J_wkday;
  int<lower=1, upper=J_station> jj_station[N];
  int<lower=1, upper=J_wkday> jj_wkday[N];
  vector[N] temp;
  vector[N] rain;
  vector[N] log_count;
  real<lower=0> scale_wip;
  real<lower=0> scale_sig;
  int<lower = 0, upper = 1> run_estimation;
}
transformed data {
  int<lower=1> K = 3;
}
parameters {
  // z refers to standardized version of the variable
  real a_z;
  real mu_a_station;
  real<lower=0> sig_a_station;
  vector[J_station] a_station_z;
  real mu_a_wkday;
  real mu_b_temp;
  real mu_b_rain;
  real<lower=0> sig_a_wkday;
  real<lower=0> sig_b_temp;
  real<lower=0> sig_b_rain;
  corr_matrix[K] Omega_wkday;
  vector[K] b_wkday[J_wkday];
  real<lower=0> sig_log_count;
}
transformed parameters {
  vector[J_station] a_station;
  real a;

  a_station = mu_a_station + sig_a_station * a_station_z;
  a = mean(log_count) + 2*sd(log_count) * a_z;
}
model {
  mu_a_station ~ normal(0, 0.01);
  sig_a_station ~ cauchy(0, scale_sig);
  a_station_z ~ normal(0, 1);

  mu_a_wkday ~ normal(0, 0.01);
  mu_b_temp ~ normal(0, scale_wip);
  mu_b_rain ~ normal(0, scale_wip);
  sig_a_wkday ~ cauchy(0, scale_sig);
  sig_b_temp ~ cauchy(0, scale_sig);
  sig_b_rain ~ cauchy(0, scale_sig);
  Omega_wkday ~ lkj_corr(2);
  {
    vector[K] mu = [mu_a_wkday, mu_b_temp, mu_b_rain]';
    vector[K] tau = [sig_a_wkday, sig_b_temp, sig_b_rain]';
    b_wkday ~ multi_normal(mu, quad_form_diag(Omega_wkday, tau));
  }

  a_z ~ normal(0, 1);
  sig_log_count ~ cauchy(0, scale_sig);
  if (run_estimation == 1) {
    vector[N] wk_i;
    for (n in 1:N)
      wk_i[n] = [1, temp[n], rain[n]] * b_wkday[jj_wkday[n]];

    log_count ~ normal(a + a_station[jj_station] + wk_i, sig_log_count);
  }
}
generated quantities {
  vector[N] y_sim;
  {
    vector[N] mu_i;
    for (n in 1:N)
      mu_i[n] = a_station[jj_station[n]] +
        [1, temp[n], rain[n]] * b_wkday[jj_wkday[n]];

    for (n in 1:N)
      y_sim[n] = normal_rng(a +  mu_i[n], sig_log_count);
  }
}
