#' Stock returns and financial conditions.
#'
#' A dataset containing the S&P 500 stock returns and the NFCI
#'
#' @format A data frame with 11,306 rows and 5 variables:
#' \describe{
#'   \item{date}{date}
#'   \item{return}{daily S&P 500 log returns times 100}
#'   \item{rv}{5-minute realized variances}
#'   \item{week}{a dummy for each year/week combination}
#'   \item{nfci}{National Financial Conditions Index}
#' }
#' @source \url{https://github.com/onnokleen/mfGARCH/}
#' @source \url{https://finance.yahoo.com/}
#' @source \url{https://fred.stlouisfed.org/series/NFCI}
#' @source \url{https://realized.oxford-man.ox.ac.uk}
"df_financial"

#' Mixed-frequency data set.
#'
#' A dataset containing the S&P 500 stock returns, realized variances and macroeconomic variables
#'
#' @format A data frame with 11,938 rows and 11 variables:
#' \describe{
#'   \item{date}{date}
#'   \item{return}{daily S&P 500 log returns times 100}
#'   \item{open_close}{open-close returns}
#'   \item{rv}{5-minute realized variances}
#'   \item{vix}{Cboe VIX}
#'   \item{year_week}{a dummy for each year/week combination}
#'   \item{dhousing}{changes in housing starts}
#'   \item{dindpro}{changes in industrial production}
#'   \item{nai}{NAI}
#'   \item{nfci}{National Financial Conditions Index}
#'   \item{year_month}{a dummy for each year/month combination}
#' }
#' @source \url{https://github.com/onnokleen/mfGARCH/}
#' @source \url{https://finance.yahoo.com/}
#' @source \url{https://fred.stlouisfed.org}
#' @source \url{https://realized.oxford-man.ox.ac.uk}
"df_mfgarch"
