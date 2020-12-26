#' @title 200-Meter Olympic Records from 1900 to 2020
#' @description A data set containing the series of record times and record
#'   values of the Olympic 200-meter, from 1900 to 2020.
#'   The variables are the following:
#'
#' \itemize{
#'   \item year  : Year of the record time
#'   \item time  : Record time
#'   \item value : Record value
#' }
#'
#' @note In this dataset, the interest lies in the lower records. The Olympics
#'   are held every 4 years, but not on all those occasions it has been held, 
#'   so those moments do not count as time.
#' @docType data
#' @keywords datasets
#' @name Olympic_records_200m
#' @source \href{https://www.olympic.org}{https://www.olympic.org}
#' @usage data(Olympic_records_200m)
#' @seealso \code{\link{series_record}}
#' @format A data frame with 12 rows and 3 variables.
"Olympic_records_200m"
