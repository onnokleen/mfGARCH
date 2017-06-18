devtools::load_all()
library(mfGARCH)
library(dplyr)

simulate_mfgarch(n.days = 1000, mu = 0, alpha = 0.02, beta = 0.9, gamma = 0, m = 0, theta = 1, w2 = 2, K = 10, psi = 0.4, sigma.psi = 0.4)


mgarch_52 <- fit_mfgarch(data = filter(df_financial, date >="1973-01-01", is.na(nfci) == FALSE),
            y = "return",
            x = "nfci",
            low.freq = "year_week",
            var.ratio.freq = "year_week",
            K = 52)
profvis::profvis(
  mgarch_0 <- fit_mfgarch(data = filter(df_demo, date >="1974-01-01"),
                          y = "return",
                          K = 0)
)


mgarch_1 <- fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
                        y = "return",
                        x = "NFCI",
                        low.freq = "week_id",
                        K = 1)


fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low.freq = "week_id",
            K = 1)

profvis::profvis(
  fit_mfgarch(data = filter(df_demo, date >="1973-01-01", is.na(nfci) == FALSE),
              y = "return",
              x = "nfci",
              low.freq = "year_week",
              K = 52)
)


profvis::profvis(
  fit_mfgarch(data = filter(df_demo, date >="1990-01-01", is.na(vix) == FALSE),
              y = "return",
              x = "vix",
              K = 3)
)

saveRDS(x, file = "test.rds")

x_test <- readRDS(file = "test.rds")

x == x_test

x_2 <- fit_mfgarch(df = dplyr::filter(df_mf_financial, Date >="1974-01-01"),
                 y = "return",
                 x = "NFCI",
                 low_freq = "week_id",
                 K = 0)
