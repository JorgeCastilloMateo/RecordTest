#' @title Extracts a subset of uncorrelated vectors
#' @importFrom stats cor.test
#' @description Given a set of \eqn{M} vectors, this function extracts 
#'   a subset of them which are uncorrelated.
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the record inference tools in this package.
#'
#'   Given a set of \eqn{M} vectors, which are the columns of matrix 
#'   \code{XM_T}, this function extracts the biggest subset of uncorrelated 
#'   vectors (columns),  using the following procedure: starting from column 
#'   \code{m}, the test \code{\link{cor.test}} is applied to study the 
#'   correlation between columns \code{m} and 
#'   \eqn{\code{m} + 1, \code{m} + 2, \ldots} an so on up to find a column 
#'   \eqn{\code{m} + k} which is not significantly correlated with column 
#'   \code{m}. Then, the process is repeated starting at column 
#'   \eqn{\code{m} + k}.
#' @param XM_T A numeric matrix (or data frame) where the uncorrelated vectors 
#'   are extracted from.
#' @param m Integer value giving the starting column.
#' @param alpha Numeric value in \eqn{(0,1)}. It gives the significance level 
#'   of the correlation test where alternative hypothesis is that the true 
#'   correlation is not equal to 0.
#' @return A vector with the index of the uncorrelated columns in the matrix.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_rev}}, 
#'   \code{\link{series_split}}, \code{\link{series_untie}}
#' @examples
#' ZM_T <- series_split(TX_Zaragoza$TX)
#' 
#' series_uncor(ZM_T)
#' 
#' @export series_uncor
series_uncor <- function(XM_T, m = 1, alpha = 0.05) {
  
  Mcols <- NCOL(XM_T)
  sep   <- c()
  while (m <= Mcols - 1) {
    sep <- c(sep, m)
    separation <- 0
    pv <- 0
    while (pv < alpha && m + separation <= Mcols - 1) {
      separation <- separation + 1
      pv <- stats::cor.test(XM_T[, m], XM_T[, m + separation])$p.value
    }
    m <- m + separation
  }
  
  # avoid correlation between the last and the first columns
  a  <- c(XM_T[, 1], NA)
  b  <- c(NA, XM_T[, sep[length(sep)]])
  pv <- stats::cor.test(a, b)$p.value
  while (pv < alpha) {
    sep <- sep[-length(sep)]
    b   <- c(NA, XM_T[, sep[length(sep)]])
    pv  <- stats::cor.test(a, b)$p.value
  }
  
  return(sep)
}