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

getSymbols.yahoo("^GSPC", env = .GlobalEnv, from = "1971-01-01")
df_returns <- as_data_frame(GSPC)
colnames(df_returns) <- c("open", "high", "low", "close", "volume", "adj_close")
df_returns$date <- as.Date(rownames(df_returns))

df_returns <-
  df_returns %>%
  dplyr::select(date, adj_close) %>%
  tbl_df()

# Download realized measures of volatility -----------------------
download.file("https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices.zip",
              destfile = "data-raw/OxfordManRealizedVolatilityIndices.zip")

df_realized <-
  read_csv("data-raw/OxfordManRealizedVolatilityIndices.zip") %>%
  filter(Symbol == ".SPX") %>%
  dplyr::select(X1 , rv5) %>%
  na.omit %>%
  mutate(date = as.Date(X1)) %>%
  rename(rv = rv5) %>%
  mutate(rv = rv * 10000) %>%
  select(date, rv)

# Join daily data -----------------------
df_daily <-
  df_returns %>%
  left_join(df_realized, by = "date") %>%
  filter(date < "1990-01-01" | is.na(rv) == FALSE) %>%
  arrange(date) %>%
  mutate(return = (log(adj_close) - log(lag(adj_close))) * 100) %>%
  mutate(rv_sq_smooth = rollapplyr(return^2, width = 22, FUN = sum, fill = NA)) %>%
  ungroup() %>%
  mutate(year = year(date), month = month(date), day = day(date)) %>%
  mutate(week = floor_date(date, unit = "weeks", week_start = 5))


df_nfci <-
  alfred::get_fred_series("NFCI", "nfci") %>%
  mutate(week = ceiling_date(date, unit = "weeks", week_start = 5)) %>%
  select(-date)

df_mf_financial <-
  df_daily %>%
  left_join(df_nfci, by = c("week")) %>%
  left_join(df_vix, by = "date") %>%
  mutate(year_month = floor_date(date, unit = "months")) %>%
  group_by(year_month) %>%
  mutate(rv_month = sum(return^2)) %>%
  ungroup() %>%
  mutate(quarter = quarter(date)) %>%
  mutate(year_quarter = floor_date(date, unit = "quarter")) %>%
  group_by(year_quarter) %>%
  mutate(rv_quarter = sum(return^2)) %>%
  ungroup() %>%
  filter(is.na(nfci) == FALSE) %>%
  select(date, return, rv, week, nfci)

rm(df_vix, df_nfci, df_daily, df_returns, df_realized, GSPC)


write_csv(df_mf_financial, "data-raw/df_mf_financial.csv")
devtools::use_data(df_mf_financial, overwrite = TRUE)

