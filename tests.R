load_all()

x <- fit_mfgarch(df = dplyr::filter(df_mf_financial, Date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low_freq = "week_id",
            K = 52)


x_2 <- fit_mfgarch(df = dplyr::filter(df_mf_financial, Date >="1974-01-01"),
                 y = "return",
                 x = "NFCI",
                 low_freq = "week_id",
                 K = 0)
