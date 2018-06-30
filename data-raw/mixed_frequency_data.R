# Mixed-frequency financial data frame -------------------------
rm(list= ls()[])
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(readr)
library(lubridate)
library(alfred)
library(quantmod)

# Download and join VIX data -----------------------

download.file("http://www.cboe.com/publish/scheduledtask/mktdata/datahouse/vixarchive.xls",
              destfile = "data-raw/VIX_1990-2003.xls")
df_vix_1990_2003 <-
  readxl::read_excel("data-raw/VIX_1990-2003.xls", skip = 1, na = "n/a") %>%
  mutate(date = as.Date(Date, "%Y-%m-%d")) %>%
  select(-Date)

df_vix_2004_current <-
  read_csv("http://www.cboe.com/publish/scheduledtask/mktdata/datahouse/vixcurrent.csv",
           skip = 1,
           col_types = "cdddd") %>%
  mutate(date = as.Date(Date, "%m/%d/%Y")) %>%
  select(-Date)

df_vix <-
  full_join(df_vix_1990_2003, df_vix_2004_current,
            by = c("VIX Open", "VIX High", "VIX Low", "VIX Close", "date"))
colnames(df_vix) <- c("VIX_open", "VIX_high", "VIX_low", "vix", "date")
df_vix <- select(df_vix, date, vix)
rm(df_vix_1990_2003, df_vix_2004_current)


# Download return data from yahoo finance -----------------------

getSymbols("^GSPC", from = "1970-01-01")
df_returns <-
  data_frame(date = index(GSPC),
             open = as.numeric(GSPC$GSPC.Open[, 1]),
             close = as.numeric(GSPC$GSPC.Close[, 1])) %>%
  arrange(date) %>%
  mutate(return = (log(close) - log(lag(close))) * 100,
         open_close = 100 * (log(close) - log(open))) %>%
  select(-open, -close)
rm(GSPC)

# Download realized measures of volatility -----------------------
download.file("https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices.zip",
              destfile = "data-raw/OxfordManRealizedVolatilityIndices.zip")
system("unzip -o data-raw/OxfordManRealizedVolatilityIndices.zip -d data-raw/")

df_realized <-
  read.csv("data-raw/oxfordmanrealizedvolatilityindices.csv") %>%
  filter(Symbol == ".SPX") %>%
  mutate(date = as.Date(substring(X, 1, 10))) %>%
  dplyr::select(date , rv5) %>%
  mutate(rv5 = as.numeric(rv5)) %>%
  rename(rv = rv5) %>%
  mutate(rv = rv * 10000) %>%
  select(date, rv)

# Join daily data -----------------------

df_daily <-
  df_returns %>%
  left_join(df_realized, by = "date") %>%
  left_join(df_vix, by = "date") %>%
  ungroup() %>%
  mutate(year_month = floor_date(date, unit = "months"),
         year_week  = floor_date(date, unit = "weeks"))
rm(df_realized, df_returns, df_vix)
saveRDS(df_daily,   file = "data-raw/df_daily.rds")

df_nai <-
  get_fred_series("CFNAI", "nai") %>%
  mutate(year_month = floor_date(date, "months")) %>%
  select(-date)
df_housing <-
  get_fred_series("HOUST", "housing") %>%
  mutate(housing = as.numeric(housing)) %>%
  mutate(dhousing = 100 * (log(housing) - lag(log(housing)))) %>%
  mutate(year_month = floor_date(date, "months")) %>%
  select(-date, -housing)
df_indpro <-
  get_fred_series("INDPRO", "indpro") %>%
  mutate(indpro = as.numeric(indpro)) %>%
  mutate(dindpro = 100 * (log(indpro) - lag(log(indpro)))) %>%
  mutate(year_month = floor_date(date, "months")) %>%
  select(-date, -indpro)
df_nfci <-
  get_fred_series("NFCI", "nfci") %>%
  mutate(year_month = floor_date(date, "months")) %>%
  mutate(year_week = floor_date(date, "weeks")) %>%
  select(-date)

df_macro <-
  list(df_housing, df_indpro, df_nai) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2, by = "year_month"), .) %>%
  arrange(year_month) %>%
  filter(year_month >= "1971-01-01")

df_mfgarch <- readRDS("data-raw/df_daily.rds") %>%
  left_join(., df_macro, by = "year_month") %>%
  select(-year_month) %>%
  left_join(., df_nfci, by = "year_week") %>%
  select(-year_month) %>%
  filter(date >= "1971-01-01") %>%
  mutate(year_month = floor_date(date, "months")) %>%
  filter(year_month <= "2018-04-01")

write_csv(df_mfgarch, "data-raw/df_mfgarch.csv")
devtools::use_data(df_mfgarch, overwrite = TRUE)

