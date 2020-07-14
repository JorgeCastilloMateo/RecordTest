#' @title Reverse elements by columns
#' @description Result from apply \code{rev} function by columns to the matrix.
#'   This allows studying the series backwards and not just forward.
#' @param XM_T A numeric matrix (or data frame).
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' series_rev(matrix(1:100, 10, 10))
#' 
#' series_rev(ZaragozaSeries)
#' 
#' @export series_rev
series_rev <- function(XM_T) apply(XM_T, 2, rev)
