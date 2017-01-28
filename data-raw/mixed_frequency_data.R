# Mixed-frequency financial data frame -------------------------
rm(list= ls()[])
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(readr)
library(lubridate)
library(Quandl)

# Download and join VIX data -----------------------

download.file("http://www.cboe.com/publish/scheduledtask/mktdata/datahouse/vixarchive.xls",
              destfile = "data-raw/VIX_1990-2003.xls")
df_vix_1990_2003 <-
  readxl::read_excel("data-raw/VIX_1990-2003.xls", skip = 1, na = "n/a") %>%
  mutate(Date = as.Date(Date, "%Y-%m-%d"))

df_vix_2004_current <-
  read_csv("http://www.cboe.com/publish/scheduledtask/mktdata/datahouse/vixcurrent.csv",
           skip = 1,
           col_types = "cdddd") %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y"))

df_vix <-
  full_join(df_vix_1990_2003, df_vix_2004_current,
            by = c("Date", "VIX Open", "VIX High", "VIX Low", "VIX Close"))
colnames(df_vix) <- c("Date", "VIX_open", "VIX_high", "VIX_low", "vix")
df_vix <- select(df_vix, Date, vix)
rm(df_vix_1990_2003, df_vix_2004_current)


# Download return data from yahoo finance -----------------------
df_returns <- Quandl("YAHOO/INDEX_GSPC", start_date = "1969-11-30")
colnames(df_returns) <- c("Date", "Open", "High", "Low", "Close", "volume", "adj_close")
df_returns <-
  df_returns %>%
  dplyr::select(Date, adj_close) %>%
  arrange(Date) %>%
  mutate(return = (log(adj_close) - log(lag(adj_close))) * 100) %>%
  mutate(rv_sq_smooth = rollapplyr(return^2, width = 22, FUN = sum, fill = NA)) %>%
  filter(Date >= "1970-01-01") %>%
  tbl_df()

# Download realized measures of volatility -----------------------
download.file("http://realized.oxford-man.ox.ac.uk/media/1366/oxfordmanrealizedvolatilityindices.zip",
              destfile = "data-raw/OxfordManRealizedVolatilityIndices.zip")

df_realized <-
  read_csv("data-raw/OxfordManRealizedVolatilityIndices.zip", skip = 2) %>%
  dplyr::select(DateID, SPX2.rv, SPX2.closeprice) %>%
  na.omit %>%
  mutate(Date = as.Date(paste(substr(DateID, start = 1, stop = 4), substr(DateID, start = 5, stop = 6), substr(DateID, start = 7, stop = 8), sep = "-"))) %>%
  rename(real_vol = SPX2.rv) %>%
  mutate(real_vol = real_vol * 10000) %>%
  select(Date, real_vol)

# Join daily data -----------------------
df_daily <-
  df_returns %>%
  left_join(df_realized, by = "Date") %>%
  ungroup() %>%
  mutate(year = year(Date), month = month(Date), day = day(Date))


df_nfci <- readxl::read_excel("data-raw/NFCI_1.xls", sheet = 2) %>%
  full_join(readxl::read_excel("data-raw/NFCI_4.xls", sheet = 2), by = c("observation_date", "NFCI")) %>%
  rename(Date = observation_date) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(year = year(Date), month = month(Date)) %>%
  select(-realtime_start_date) %>%
  mutate(week_id = seq(1:dim(.)[1]))



df_mf_financial <- df_daily %>%
  left_join(df_nfci, by = c("year", "month", "Date")) %>%
  left_join(df_vix, by = "Date") %>%
  mutate(year_month = as.numeric(paste0(year, sprintf("%02d", .$month)))) %>%
  group_by(year_month) %>%
  mutate(rv.month = sum(return^2)) %>%
  ungroup() %>%
  mutate(quarter = quarter(Date)) %>%
  mutate(year_quarter = as.numeric(paste0(year, sprintf("%02d", .$quarter)))) %>%
  group_by(year_quarter) %>%
  mutate(rv_quarter = sum(return^2)) %>%
  filter(is.na(real_vol) == FALSE || Date < "2000-01-01") %>%
  ungroup()

rm(df_vix, df_nfci, df_daily, df_returns, df_realized)

last_NFCI <- NA
last_week_id <- NA
# Fill in weekly data
for (i in c(1:dim(df_mf_financial)[1])) {
  if (is.na(df_mf_financial$NFCI[i]) == TRUE) {
    df_mf_financial$NFCI[i] <- last_NFCI
    df_mf_financial$week_id[i] <- last_week_id
  } else{
    last_NFCI <- df_mf_financial$NFCI[i]
    last_week_id <- df_mf_financial$week_id[i]
  }
}
rm(i, last_NFCI, last_week_id)

write_csv(df_mf_financial, "data-raw/df_mf_financial.csv")
devtools::use_data(df_mf_financial, overwrite = TRUE)
