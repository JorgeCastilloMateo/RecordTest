#' @title The Classical Record Model
#' @importFrom stats rnorm
#' @description Random generation for the classical record model, i.e.,
#'   sequences of independent and identically distributed (IID) continuous
#'   random variables (RVs).  
#' @param Trows,Mcols Integers indicating the number of rows and columns of the
#'   returned matrix, i.e., the length and number of series for the record 
#'   analysis.
#' @param rdist A function that simulates continuous random variables, 
#'   e.g., \code{\link{runif}} (fastest in \code{stats} package), 
#'   \code{\link{rnorm}} or \code{\link{rexp}}.
#' @param ... Further arguments to introduce in the \code{rdist} function.
#' @return A matrix of draws of IID continuous RVs with common distribution 
#'   \code{rdist}.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L.record}}, \code{\link{S.record}}, 
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}, 
#'   \code{\link{p.record}}, \code{\link{records}}
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. 
#' Wiley Series in Probability and Statistics. Wiley, New York.
#' @examples
#' # By default, draw a sample of 100 series of length 50 
#' # with observations coming from a standard normal distribution 
#' X <- rcrm()
#' # Compute its record indicators
#' I <- I.record(X)
#' # Implement some tests
#' N.test(X, distribution = "poisson-binomial")
#' foster.test(X, weights = function(t) t-1, statistic = "D")
#' 
#' @export rcrm
rcrm <- function(Trows = 50, Mcols = 100, rdist = stats::rnorm, ...) {
  
  return(matrix(rdist(Trows * Mcols, ...), nrow = Trows, ncol = Mcols))
}
