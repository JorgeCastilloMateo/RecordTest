#' @title 200-Meter Olympic Records from 1900 to 2020
#' @description A data set containing the record times and record values of the
#'   200-meter competition at the Olympic games, from 1900 to 2020.
#'   The variables are the following:
#'
#' \itemize{
#'   \item year  : Year of the record time
#'   \item time  : Record time
#'   \item value : Record value in seconds
#' }
#'
#' @note In this data set, the interest lies in the lower records. Although 
#'   the Olympic Games are held every 4 years, not all of these occasions have
#'   been held, so only the games that have taken place are considered in the 
#'   definition of time.
#' @docType data
#' @keywords datasets
#' @name Olympic_records_200m
#' @source \href{https://www.olympic.org}{https://www.olympic.org}
#' @usage data(Olympic_records_200m)
#' @seealso \code{\link{series_record}}
#' @format A data frame with 12 rows and 3 variables.
"Olympic_records_200m"
