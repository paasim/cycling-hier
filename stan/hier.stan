data {
  int<lower=0> N;
  int<lower=1> J_station;
  int<lower=1> J_wkday;
  int<lower=1, upper=J_station> jj_station[N];
  int<lower=1, upper=J_wkday> jj_wkday[N];
  vector[N] temp;
  vector[N] rain;
  vector<lower=1>[N] count;
  int<lower=0, upper=1> run_estimation;
}
transformed data {
  vector<lower=0, upper=1>[N] is_rain;
  vector<lower=0>[N] y = log(count);
  for (n in 1:N)
    is_rain[n] = (rain[n] > 0) ? 1 : 0;
}
parameters {
  // z refers to standardized version of the variable
  vector[J_station] a_station_z;
  vector[J_wkday] a_wkday_z;
  vector[J_wkday] a_rain_z;
  vector[J_wkday] b_temp_z;
  real mu_a_station;
  real mu_a_wkday;
  real mu_a_rain;
  real mu_b_temp;
  real<lower=0> sig_a_station;
  real<lower=0> sig_a_wkday;
  real<lower=0> sig_a_rain;
  real<lower=0> sig_b_temp;

  real a_z;
  real<lower=0> sig_y;
}
transformed parameters {
  vector[J_station] a_station = mu_a_station + sig_a_station * a_station_z;
  vector[J_wkday] a_wkday = mu_a_wkday + sig_a_wkday * a_wkday_z;
  vector[J_wkday] a_rain = mu_a_rain + sig_a_rain * a_rain_z;
  vector[J_wkday] b_temp  = mu_b_temp + sig_b_temp * b_temp_z;
  real a = mean(y) + 2*sd(y) * a_z;
}
model {
  a_station_z ~ normal(0, 1);
  a_wkday_z ~ normal(0, 1);
  a_rain_z ~ normal(0, 1);
  b_temp_z ~ normal(0, 1);

  mu_a_station ~ normal(0, 0.01);
  mu_a_wkday ~ normal(0, 0.01);
  mu_a_rain ~ normal(0, 1);
  mu_b_temp ~ normal(0, 1);
  sig_a_station ~ cauchy(0, 2.5);
  sig_a_wkday ~ cauchy(0, 2.5);
  sig_a_rain ~ normal(0, 0.5);
  sig_b_temp ~ normal(0, 0.5);

  a_z ~ normal(0, 1);
  sig_y ~ cauchy(0, 2.5);
  if (run_estimation == 1) {
    vector[N] mu_y = a + a_station[jj_station] + a_wkday[jj_wkday] +
                      a_rain[jj_wkday] .* is_rain + b_temp[jj_wkday] .* temp;

    y ~ normal(mu_y, sig_y);
  }
}
generated quantities {
  vector[N] count_sim;
  {
    vector[N] sim;
    vector[N] mu_y = a + a_station[jj_station] + a_wkday[jj_wkday] +
                      a_rain[jj_wkday] .* is_rain + b_temp[jj_wkday] .* temp;
    for (n in 1:N)
      sim[n] = normal_rng(0, 1);

    count_sim = exp(mu_y + sim * sig_y);
  }
}
