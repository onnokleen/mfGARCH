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
#' @source \url{https://de.finance.yahoo.com/}
#' @source \url{https://fred.stlouisfed.org/series/NFCI}
#' @source \url{https://realized.oxford-man.ox.ac.uk}
"df_financial"
