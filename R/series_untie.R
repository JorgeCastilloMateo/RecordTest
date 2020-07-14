#' @title Avoid record ties
#' @description Avoids record ties when observations have been rounded.
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param a,b Lower and upper limits of the uniform distribution.
#' @param seed Allows the user to specify the seed of the random sampling. 
#'   If \code{NULL}, then no seed is specified.
#' @details This function is used in the data preparation (or pre-processing). 
#'  Its use is recommended but not necessary to apply the record inference 
#'  tools in this package.
#'  
#'  The record theory on which the hypothesis tests are based assumes 
#'  that the random variables are continuous, proving that the probability 
#'  that two observations take the same value is zero. Most of the data 
#'  collected is rounded, giving a certain probability to the tie between 
#'  records, slightly skewing the results, in general obtaining a fewer 
#'  amount of records if it has been rounded a lot. 
#'  
#'  This function avoids ties by adding a sample from a uniform variable 
#'  in the interval \eqn{(\code{a}, \code{b})} to each element of \code{XM_T}.
#'  If rounding is in units it is recommended to take \code{a = -0.5, b = 0.5},
#'  while if it is in tenths it is recommended \code{a = -0.05, b = 0.05}
#'  and so on.
#'  
#'  If the function is used correctly (appropiate \code{a,b} values) the records
#'  in the original series will also be in the new series.
#'
#' @return A matrix equal to \code{XM_T} whose elements have been added a 
#'   sample from a uniform variable.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_rev}}
#'   \code{\link{series_split}}, \code{\link{series_uncor}} 
#' @examples
#' set.seed(23)
#' series_untie(matrix(round(rnorm(100), digits = 1), 10, 10), a = -0.05, b = 0.05)
#' 
#' series_untie(ZaragozaSeries, seed = 23)
#' 
#' @export series_untie
series_untie <- function(XM_T, a = -0.5, b = 0.5, seed = NULL) {
  
  if (!is.null(seed)) {
    old <- .Random.seed
    on.exit( { .Random.seed <<- old } )
    set.seed(seed)
  }
  
  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  
  XM_T <- XM_T + matrix(runif(Trows * Mcols, min = a, max = b), nrow = Trows, ncol = Mcols)
  
  return(XM_T)
}
