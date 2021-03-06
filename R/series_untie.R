#' @title Avoid Record Ties
#' @importFrom stats runif rnorm
#' @description Avoids record ties when observations have been rounded.
#' @param X A numeric vector, matrix (or data frame).
#' @details This function is used in the data preparation (or pre-processing) 
#'  often required to apply the exploratory and inference tools based on 
#'  theory of records within this package.
#'  
#'  The theory of records on which the hypothesis tests are based assumes 
#'  that the random variables are continuous, proving that the probability 
#'  that two observations take the same value is zero. Most of the data 
#'  collected is rounded, giving a certain probability to the tie between 
#'  records, thereby reducing the number of new records (see, e.g., Wergen 
#'  et al. 2012). 
#'  
#'  This function avoids ties by adding a sample from a uniform random variable 
#'  to each element of \code{X}. The function computes the maximum number of 
#'  decimal digits (let it be \eqn{n}) for any element in \code{X}. Then the
#'  uniform variable is sampled in the interval 
#'  \eqn{(-5 \times 10^{-(n+1)}, 5 \times 10^{-(n+1)})}, so the records in the
#'  original (rounded) series will also be in the new series.
#'
#' @return A matrix equal to \code{X} whose elements have been added a 
#'   sample from a uniform variable, different for every element.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_record}}, 
#'   \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_ties}}, \code{\link{series_uncor}} 
#' @references 
#' Wergen G, Volovik D, Redner S, Krug J (2012). 
#' “Rounding Effects in Record Statistics.”
#' \emph{Physical Review Letters}, \strong{109}(16), 164102. 
#' @examples
#' set.seed(23)
#' X <- matrix(round(stats::rnorm(100), digits = 1), nrow = 10, ncol = 10)
#' series_untie(X)
#' 
#' series_untie(ZaragozaSeries)
#' 
#' @export series_untie
series_untie <- function(X) {
  
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  
  digits <- max(nchar(sub("^0+", "", sub("\\.", "", X %% 1)))) + 1
  a <- -5 * 10^-digits
  b <-  5 * 10^-digits
  
  X <- X + matrix(stats::runif(Trows * Mcols, min = a, max = b), nrow = Trows, ncol = Mcols)
  
  return(X)
}
