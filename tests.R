load_all()


simulate_mfgarch(n.days = 1000, mu = 0, alpha = 0.02, beta = 0.9, gamma = 0, m = 0, theta = 1, w2 = 2, K = 10, psi = 0.4, sigma.psi = 0.4)


mgarch_52 <- fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low.freq = "week_id",
            K = 52)

mgarch_0 <- fit_mfgarch(df = dplyr::filter(dplyr::mutate(df_mf_financial, date = Date), date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low.freq = "week_id",
            K = 0)

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
