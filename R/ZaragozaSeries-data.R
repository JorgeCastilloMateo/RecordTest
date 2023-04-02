#' @title Split and Uncorrelated Time Series \code{\link{TX_Zaragoza}}
#' 
#' @description The matrix resulting from the data preparation (or 
#'   pre-processing) of \code{\link{TX_Zaragoza}$TX}. 
#' 
#' @details 
#'   The matrix is the result from applying:
#'   \code{\link{series_uncor}(\link{series_split}(\link{TX_Zaragoza}$TX))}.
#'   
#'   The data matrix corresponds to the 70 years with observations in 
#'   \code{\link{TX_Zaragoza}$TX} and to the 76 days in the year where adjacent
#'   daily maximum temperature subseries are uncorrelated. By coincidence, 
#'   none of the subseries 4, 90 or 278 with missing observations is kept 
#'   within the 76 uncorrelated days.
#'
#' @docType data
#' @keywords datasets
#' @name ZaragozaSeries
#' @seealso \code{\link{TX_Zaragoza}}
#' @usage data("ZaragozaSeries")
#' @format A matrix with 70 rows and 76 columns.
"ZaragozaSeries"
