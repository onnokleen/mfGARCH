load_all()

fit_mfgarch(df = df_mf_financial,
            y = "return",
            x = "NFCI",
            low_freq = "week_id",
            K = 12)
