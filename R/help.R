#'@section RequestLogs:
#'We don't study readers enough, and we finally have the tools to do that. WMUtils contains
#'several functions centred on the RequestLogs. \code{\link{hive_query}} allows you to query
#'the unsampled logs, while \code{\link{sampled_logs}} allows you to retrieve the 1:1000 sampled
#'ones. For both data types, \code{\link{log_strptime}} turns the timestamp format used into
#'POSIXlt timestamps.
#'
#'@section MySQL:
#'If you study editors, our MySQL databases are where all the data lives. \code{\link{mysql_query}}
#'allows you to query a single database on \code{analytics-store.eqiad.wmnet}, while
#'\code{\link{global_query}} allows you to run over multiple databases. Either way,
#'\code{\link{mw_strptime}} turns the timestamp format used in our DB into POSIXlt timestamps.
#'And once you're done processing, use \code{\link{mysql_write}} to stream the results up to
#'the databases again. Need to update previously written rows? No problem! \code{\link{mysql_delete}}
#'is the function for you.


#' @useDynLib mfGARCH
#' @importFrom Rcpp sourceCpp
NULL
