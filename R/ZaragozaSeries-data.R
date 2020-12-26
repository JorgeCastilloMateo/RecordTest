#' @title Splitted and Uncorrelated Time Series \code{\link{TX_Zaragoza}}
#' @description The matrix resulting from the pre-processing of 
#'   \code{\link{TX_Zaragoza}$TX}. 
#' 
#' @details 
#'   The matrix is the result from applying:
#'   \code{\link{series_uncor}(\link{series_split}(\link{TX_Zaragoza}$TX))}.
#'
#' @docType data
#' @keywords datasets
#' @name ZaragozaSeries
#' @usage data(ZaragozaSeries)
#' @seealso \code{\link{TX_Zaragoza}}
#' @format A matrix with 66 rows and 76 columns.
"ZaragozaSeries"
