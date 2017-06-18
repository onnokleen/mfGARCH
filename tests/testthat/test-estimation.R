test_that("Error testing", {
  expect_warning(
    mgarch_0 <- fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01"),
                            y = "return",
                            K = 0)$par
  )

  expect_warning(
    mgarch_0 <- fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01"),
                            y = "return",
                            gamma = FALSE,
                            K = 0)$par
  )

  expect_equal(
    fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01", is.na(nfci) == FALSE),
                y = "return",
                x = "nfci",
                low.freq = "year_week",
                K = 52)$variance.ratio,
    13.1500626717928
  )

  expect_equal(
    fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01", is.na(nfci) == FALSE),
                y = "return",
                x = "nfci",
                low.freq = "year_week",
                gamma = FALSE,
                K = 52)$variance.ratio,
    10.4571948652477
  )

  expect_error( # K should not be smaller than 0 and should be an integer
    fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01"),
                y = "return",
                x = "nfci",
                K = -10)
  )
  expect_error( # K should not be smaller than 0 and should be an integer
    fit_mfgarch(data = dplyr::filter(df_financial, date >="1974-01-01"),
                y = "return",
                x = "nfci",
                K = 0.22)
  )
  expect_error( # nfci includes NAs
    fit_mfgarch(data = df_financial,
                y = "return",
                x = "nfci",
                K = 52)
  )
  expect_error(
    fit_mfgarch(data = df_financial,
                y = "return",
                x = "nfci",
                K = 52)
  )
  # expect_equal(
  #   get_alfred_series("INDPRO", "test",
  #                     observation_start = "2013-03-01", observation_end = "2013-03-01",
  #                     realtime_start = "2015-02-02", realtime_end = "2015-02-02"),
  #   data.frame(
  #     date = as.Date("2013-03-01"),
  #     realtime_period  = as.Date("2015-02-02"),
  #     test = 99.488)
  # )
  # expect_equal(
  #   dplyr::filter(get_alfred_series("INDPRO", "test",
  #                                   realtime_start = "2015-02-02", realtime_end = "2015-02-02"),
  #                 date == "2013-03-01"),
  #   data.frame(
  #     date = as.Date("2013-03-01"),
  #     realtime_period  = as.Date("2015-02-02"),
  #     test = 99.488)
  # )
  # expect_equal(
  #   dplyr::filter(get_alfred_series("INDPRO", "test",
  #                                   observation_start = "2009-03-01", observation_end = "2009-03-01"),
  #                 realtime_period == "2015-02-18"),
  #   data.frame(
  #     date = as.Date("2009-03-01"),
  #     realtime_period = as.Date("2015-02-18"),
  #     test = 85.6157)
  # )
  # expect_equal(
  #   get_fred_series("INDPRO", "test",
  #                   observation_start = "2009-03-01", observation_end = "2009-03-01"),
  #   data.frame(
  #     date = as.Date("2009-03-01"),
  #     test = 89.1913)
  # )
  # expect_equal(
  #   dplyr::filter(get_fred_series("INDPRO", "test"),
  #                 date == "2009-03-01"),
  #   data.frame(
  #     date = as.Date("2009-03-01"),
  #     test = 89.1913)
  # )
  # expect_error(
  #   get_fred_series(1231232, "test",
  #                   observation_start = "2009-03-01", observation_end = "2009-03-01")
  # )
  # expect_error(
  #   get_alfred_series(1231232, "test",
  #                     observation_start = "2009-03-01", observation_end = "2009-03-01")
  # )
})
