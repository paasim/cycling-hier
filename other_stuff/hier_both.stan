data {
  int<lower=0> N;
  int<lower=1> J_station;
  int<lower=1> J_wkday;
  int<lower=1, upper=J_station> jj_station[N];
  int<lower=1, upper=J_wkday> jj_wkday[N];
  vector[N] temp;
  vector[N] rain;
  vector[N] log_count;
  real scale_wip;
  real scale_sig;
  int<lower = 0, upper = 1> run_estimation;
}
transformed data {
  vector<lower=0, upper=1>[N] is_rain;
  for (n in 1:N)
    is_rain[n] = (rain[n] == 0) ? 0 : 1;
}
parameters {
  // z refers to standardized version of the variable
  real a_z;
  vector[J_station] a_station_z;
  real mu_a_station;
  real<lower=0> sig_a_station;
  real mu_a_wkday;
  real mu_b_temp;
  real mu_a_rain;
  real mu_b_rain;
  real<lower=0> sig_a_wkday;
  real<lower=0> sig_b_temp;
  real<lower=0> sig_a_rain;
  real<lower=0> sig_b_rain;
  vector[J_wkday] a_wkday_z;
  vector[J_wkday] b_temp_z;
  vector[J_wkday] a_rain_z;
  vector[J_wkday] b_rain_z;
  real<lower=0> sig_log_count;
}
transformed parameters {
  vector[J_station] a_station;
  vector[J_wkday] a_wkday;
  vector[J_wkday] b_temp;
  vector[J_wkday] a_rain;
  vector[J_wkday] b_rain;
  real a;

  a_station = mu_a_station + sig_a_station * a_station_z;

  a_wkday = mu_a_wkday + sig_a_wkday * a_wkday_z;
  b_temp = mu_b_temp + sig_b_temp * b_temp_z;
  a_rain = mu_a_rain + sig_a_rain * a_rain_z;
  b_rain = mu_b_rain + sig_b_rain * b_rain_z;

  a = mean(log_count) + 2*sd(log_count) * a_z;
}
model {
  mu_a_station ~ normal(0, 0.01);
  sig_a_station ~ cauchy(0, scale_sig);
  a_station_z ~ normal(0, 1);

  mu_a_wkday ~ normal(0, 0.01);
  mu_b_temp ~ normal(0, scale_wip);
  mu_a_rain ~ normal(0, scale_wip);
  mu_b_rain ~ normal(0, scale_wip);
  sig_a_wkday ~ cauchy(0, scale_sig);
  sig_b_temp ~ cauchy(0, scale_sig);
  sig_a_rain ~ cauchy(0, scale_sig);
  sig_b_rain ~ cauchy(0, scale_sig);
  a_wkday_z ~ normal(0, 1);
  b_temp_z ~ normal(0, 1);
  a_rain_z ~ normal(0, 1);
  b_rain_z ~ normal(0, 1);

  a_z ~ normal(0, 1);
  sig_log_count ~ cauchy(0, scale_wip);
  if (run_estimation == 1) {
    vector[N] mu_i = a + a_station[jj_station] + a_wkday[jj_wkday] +
                      a_rain[jj_wkday] .* is_rain +
                      b_temp[jj_wkday] .* temp + b_rain[jj_wkday] .* rain;

    log_count ~ normal(mu_i, sig_log_count);
  }
}
generated quantities {
  vector[N] y_sim;
  {
    vector[N] mu_i = a + a_station[jj_station] + a_wkday[jj_wkday] +
                      a_rain[jj_wkday] .* is_rain +
                      b_temp[jj_wkday] .* temp + b_rain[jj_wkday] .* rain;
    for (n in 1:N)
      y_sim[n] = normal_rng(mu_i[n], sig_log_count);
  }
}
