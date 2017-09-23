# simplefmi is available at github: devtools::install_github("paasim/simplefmi")
# also, to use the fmi weather api, an api key is required, see
# http://en.ilmatieteenlaitos.fi/open-data-manual-fmi-wfs-services
library(tidyverse)
library(forcats)
library(stringr)
library(feather)
library(lubridate)
library(simplefmi)


# Read the cycling data
data_url <- "http://www.hel.fi/hel2/tietokeskus/data/helsinki/ksv/Helsingin_pyorailijamaarat.csv"

df_raw <- read_csv2(data_url, col_names = TRUE,
                 col_types = cols("Päivämäärä" = "c", .default = "i"),
                 locale = locale(encoding = "ISO-8859-1", decimal_mark = ","))
colnames(df_raw) <- colnames(df_raw) %>% tolower() %>%
  str_replace_all(setNames(c("a", "o", "", "i ", "_"),
                           c("ä", "ö", "\\.| \\(.*\\)|(ela|hjois|a|nsi)puoli| silta", "en ", " ")))

df_2014 <- filter(df_raw, as.numeric(str_extract(paivamaara, "(?=\\D*\\d+\\D*)\\d{4}")) >= 2014)
  # exclude every observation before 2014 as the format
  # before that is not consistent with other observations

days_en <- locale()$date_names$day_ab[c(2:7,1)]
mon_fi <- locale("fi")$date_names$mon %>% str_sub(end = -6)
mon_map <- function(x) str_replace_all(x, setNames(as.character(1:12), mon_fi))
varnames <- c("wkday", "day", "mon", "year", "hr")

# get a tibble with dates
df_dates <- separate(df_2014, "paivamaara", varnames, sep = " ") %>%
  mutate(date = make_date(year, mon_map(mon), day)) %>%
  select(date, huopalahti:baana) %>%
  group_by(date) %>% # get daily counts
  summarise_all(sum) %>%
  ungroup() %>%
  # combine observations that correspond to same location
  mutate(kulosaari = kulosaari_et + kulosaari_po,
         munkkiniemi = munkkiniemi_et + munkkiniemi_po,
         pitkasilta = pitkasilta_it + pitkasilta_la,
         lauttasaari = lauttasaari_et + lauttasaari_po,
         wkday = weekdays(date) %>% str_sub(end = 3) %>%
           as_factor() %>% fct_relevel(days_en)) %>%
  # drop redundant variables + hesperian puisto to get rid of zeros
  select(-matches("_\\w")) %>%
  na.omit()# drop NAs for now

# get the respective weather data:
fmi_apikey <- readLines("apik") # fmi-apikey
# station id for kaisaniemi from
# http://en.ilmatieteenlaitos.fi/observation-stations
station_id <- "100971"

df_weather <- fmi_download(fmi_apikey, min(df_dates$date), max(df_dates$date),
                           station_id, hourly = FALSE) %>%
  mutate(rain = pmax(rain, 0)) # rain = -1 => no rain.


# combine the data frames
left_join(df_dates, df_weather, "date") %>%
  # reorder the columns
  select(date, wkday, rain, temp, huopalahti:lauttasaari) %>%
  write_feather("data/df.feather")
