#' @title Subset of Uncorrelated Series
#' 
#' @importFrom stats cor.test
#' 
#' @description Given a matrix this function extracts a subset of uncorrelated
#'   columns (see Details).
#'   
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the exploratory and inference tools based on 
#'   theory of records within this package.
#'
#'   Given a matrix \code{X} considered as a set of \eqn{M^*} vectors, which
#'   are the columns of \code{X}, this function extracts a subset 
#'   of uncorrelated vectors (columns), using the following procedure: starting
#'   from column \code{m}, the test \code{test.fun} is applied to
#'   study the correlation between columns depending on argument \code{type}.
#'   
#'   If \code{type = "adjacent"}, the test is computed between \code{m}
#'   and \eqn{\code{m} + 1, \code{m} + 2, \ldots} and so on up to find a column 
#'   \eqn{\code{m} + k} which is not significantly correlated with column 
#'   \code{m}. Then, the process is repeated starting at column 
#'   \eqn{\code{m} + k}. All columns are checked. 
#'   
#'   When the first and last columns may not have a significant correlation, 
#'   where \code{m} is the first column, the parameter \code{first.last} should
#'   be \code{FALSE}. When the first and last columns could be correlated, 
#'   the function requires \code{first.last = TRUE}.
#'   
#'   If \code{type = "all"}, the procedure is similar as above but the new kept
#'   column cannot be significantly correlated with any other column already 
#'   kept, not only the previous one. So this option results in a fewer number
#'   of columns.
#'   
#' @param X A numeric matrix (or data frame) where the uncorrelated vectors 
#'   are extracted from.
#' @param test.fun A function that tests the correlation (it could also
#'   test dependence or other feature that is desired to test on the columns). 
#'   It must take as arguments two numeric vectors of the same length to apply 
#'   the test to. The alternative hypothesis between both vectors is that the 
#'   true correlation is not equal to 0 (or that they are dependent, etc). The 
#'   return value should be a list object with a component \code{p.value} with 
#'   the p-value of the test. Default is \code{\link[stats]{cor.test}}. Other 
#'   function to test tail dependence is found in the package \code{extRemes}
#'   as \code{taildep.test}.
#' @param return.value A character string indicating the return of the function,
#'  \code{"series"} for a matrix with uncorrelated columns or \code{"indexes"}
#'  for a vector with the position of the uncorrelated columns in \code{X}.
#' @param type A character string indicating the type of uncorrelation wanted
#'   between the extracted series (or columns), \code{"adjacent"} or 
#'   \code{"all"} (see Details).
#' @param first.last Logical. Indicates if the first and last columns have also
#'   to be uncorrelated (when \code{type = "adjacent"}).
#' @param m Integer value giving the starting column.
#' @param alpha Numeric value in \eqn{(0,1)}. It gives the significance level. 
#'   For \code{\link[stats]{cor.test}} the alternative hypothesis is that the 
#'   true correlation is not equal to 0.
#' @param ... Further arguments to be passed to \code{test.fun} function.
#' @return A matrix or a vector as specified by \code{return.value}.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_record}},
#'   \code{\link{series_rev}}, \code{\link{series_split}},
#'   \code{\link{series_ties}}, \code{\link{series_untie}}
#' @references 
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2), 313-330. 
#' \doi{10.1007/s00477-021-02122-w}.
#' 
#' @examples
#' # Split Zaragoza series
#' TxZ <- series_split(TX_Zaragoza$TX)
#' 
#' # Index of uncorrelated columns depending on the criteria
#' series_uncor(TxZ, return.value = "indexes", type = "adjacent")
#' series_uncor(TxZ, return.value = "indexes", type = "all")
#' 
#' # Return the set of uncorrelated vectors
#' ZaragozaSeries <- series_uncor(TxZ)
#' 
#' @export series_uncor
series_uncor <- function(X, 
                         test.fun = stats::cor.test,
                         return.value = c("series", "indexes"), 
                         type = c("adjacent", "all"),
                         first.last = TRUE,
                         m = 1, 
                         alpha = 0.05, 
                         ...) {
  
  return.value <- match.arg(return.value)
  type <- match.arg(type)
  if (!(0 < alpha & alpha < 1)) { stop("'alpha' should be in (0,1)") }
  if(!is.function(test.fun)) stop("'test.fun' should be a function")
  
  m0 <- m
  Mcols <- ncol(X)
  if (m0 != 1) { X <- X[, c(m0:Mcols, 1:(m0 - 1))] }
  sep <- m <- 1
  
  if (type == 'adjacent') {
    
    while (m < Mcols) {
      separation <- 0
      pv <- 0
      while (pv < alpha && m + separation < Mcols) {
        separation <- separation + 1
        pv <- test.fun(X[, m], X[, m + separation], ...)$p.value
      }
      if (m + separation == Mcols && pv < alpha) { break }
      m <- m + separation
      sep <- c(sep, m)
    }
    
    # avoid correlation between the first and the last columns
    if (first.last) {
      
      if (length(sep) != 1) {
        
        a  <- X[, 1]
        b  <- X[, sep[length(sep)]]
        pv <- test.fun(a, b, ...)$p.value
        
        while (pv < alpha) {
          sep <- sep[-length(sep)]
          b   <- X[, sep[length(sep)]]
          pv  <- test.fun(a, b, ...)$p.value
        }
      }
    }
    
  } else { # type == 'all'
    
    while (m < Mcols) {
      separation <- 0
      pv <- 0
      while (pv < alpha && m + separation < Mcols) {
        separation <- separation + 1
        for(i in seq_along(sep)) { pv[i] <- test.fun(X[, sep[i]], X[, m + separation], ...)$p.value }
        pv <- min(pv)
      }
      if (m + separation == Mcols && pv < alpha) { break }
      m <- m + separation
      sep <- c(sep, m)
    }
  }
  
  if (m0 != 1) {
    
    sep <- (sep + m0 - 1) %% Mcols
    sep[sep == 0] <- Mcols
    sep <- sort(sep)
  } 
  
  if (return.value == "indexes") { return(sep) }

  if (m0 != 1) { X <- X[, c((Mcols - m0 + 2):Mcols, 1:(Mcols - m0 + 1))] }
  
  return(X[, sep])
}
