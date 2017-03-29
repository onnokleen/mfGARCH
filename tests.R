load_all()

fit_mfgarch(df = dplyr::filter(df_mf_financial, Date >="1974-01-01"),
            y = "return",
            x = "NFCI",
            low_freq = "week_id",
            K = 52)

sourceCpp("src/calculate_g.cpp")
sourceCpp("src/sum_tau.cpp")
