rm(list = ls())
library(tidyverse)
library(forcats)
library(feather)
library(readr)
library(bayesplot)
library(lubridate)
library(loo)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(1)
(df_raw <- read_feather("data/df.feather"))

df <- gather(df_raw, station, count,
             -(1:5), factor_key = TRUE) %>%
  mutate(log_count = log(count))

# sample 50 random dates and
# select the corresponding observations from each station
df_tr <- filter(df, date %in% sample(df_raw$date, 160))

d <- list(N = nrow(df_tr),
           J_station = levels(df_tr$station) %>% length(),
           J_wkday = levels(df_tr$wkday) %>% length(),
           jj_station = as.integer(df_tr$station),
           jj_wkday = as.integer(df_tr$wkday),
           temp = df_tr$temp,
           rain = df_tr$rain,
           log_count = df_tr$log_count,
           scale_wip = 2,
           scale_sig = 1.5,
           run_estimation = 1)

m1 <- stan_model("stan/hier.stan")
m2 <- stan_model("stan_old/hier_coef.stan")
m3 <- stan_model("stan_old/hier_const.stan")

fit1 <- sampling(m1, data = d)
fit2 <- sampling(m2, data = d)
fit3 <- sampling(m3, data = d)

loo1 <- loo(extract(fit1)$y_sim)
loo2 <- loo(extract(fit2)$y_sim)
loo3 <- loo(extract(fit3)$y_sim)
compare(loo1, loo2, loo3)
