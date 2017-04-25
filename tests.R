load_all()
profvis(
x <- fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low.freq = "week_id",
            K = 52)
)

profvis(
  fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
              y = "return",
              x = "NFCI",
              low.freq = "week_id",
              K = 52)
)

saveRDS(x, file = "test.rds")

x_test <- readRDS(file = "test.rds")

x == x_test

x_2 <- fit_mfgarch(df = dplyr::filter(df_mf_financial, Date >="1974-01-01"),
                 y = "return",
                 x = "NFCI",
                 low_freq = "week_id",
                 K = 0)
