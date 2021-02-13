#' @title Splitted and Uncorrelated Time Series \code{\link{TX_Zaragoza}}
#' @description The matrix resulting from the data-preparation (or 
#'   pre-processing) of \code{\link{TX_Zaragoza}$TX}. 
#' 
#' @details 
#'   The matrix is the result from applying:
#'   \code{\link{series_uncor}(\link{series_split}(\link{TX_Zaragoza}$TX))}.
#'   
#'   The data matrix corresponds to the 70 years with observations in 
#'   \code{\link{TX_Zaragoza}$TX} and to the 76 days in the year where adjacent
#'   daily maximum temperature sub-series are uncorrelated. Casually, none of 
#'   the sub-series 4, 90 or 278 with imputed values is kept within the 76 
#'   uncorrelated days.
#'
#' @docType data
#' @keywords datasets
#' @name ZaragozaSeries
#' @usage data(ZaragozaSeries)
#' @seealso \code{\link{TX_Zaragoza}}
#' @format A matrix with 70 rows and 76 columns.
"ZaragozaSeries"
