#' @title Reverse Elements by Columns
#' @description Result of applying \code{\link[base]{rev}} function by columns
#'   to the matrix. This allows the study of the series backwards and not only 
#'   forward.
#' @param X A numeric vector, matrix (or data frame).
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_record}}, 
#'   \code{\link{series_split}}, \code{\link{series_ties}},
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' series_rev(matrix(1:100, 10, 10))
#' 
#' series_rev(ZaragozaSeries)
#' 
#' @export series_rev
series_rev <- function(X) { return(apply(as.matrix(X), 2, rev)) }
