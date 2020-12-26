#' @title Global Statistic for Two-Sided Tests
#' @importFrom stats runif
#' @description This function performs a more powerful generalization 
#'   of the two-sided tests in this package by means of the sum of the 
#'   statistics of upper and lower records in the forward and backward 
#'   directions to study the hypothesis of the classical record model. The
#'   tests considered are the regression test \code{\link{p.test}},
#'   the chi-square goodness-of-fit test \code{\link{chisq.test}},
#'   the likelihood-ratio test \code{\link{lr.test}} and the score test
#'   \code{\link{score.test}}.
#' @details 
#'   The statistics, say \eqn{X}, of the tests \code{\link{p.test}}, 
#'   \code{\link{chisq.test}}, \code{\link{lr.test}} or
#'   \code{\link{score.test}} for the forward upper, forward lower, backward
#'   upper and backward lower records are summed to develop a more powerful 
#'   statistic:
#'   \deqn{X^{G} = X^{(FU)} + X^{(FL)} + X^{(BU)} + X^{(BL)},}
#'   where those are the statistics for their respective type of record. 
#'   Other sums of statistics are allowed.
#'   
#'   The distribution of this global statistics is unknown, 
#'   but the p-value can be estimated with Monte Carlo simulations
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param FUN One of the functions whose statistic you want to use. One of
#'   \code{\link{chisq.test}}, \code{\link{lr.test}}, 
#'   \code{\link{p.test}} or \code{\link{score.test}}.
#' @param record Logical vector. Vector with four elements indicating if 
#'   forward upper, forward lower, backward upper and backward lower are going
#'   to be shown, respectively. Logical values or 0,1 values are accepted.
#' @param B An integer specifying the number of replicates used in the 
#'   Monte Carlo approach.
#' @param ... Further arguments in the \code{FUN} function.
#' @return A list of class \code{"htest"}  with the following elements:
#'   \item{statistic}{Value of the  statistic.}
#'   \item{p.value}{Simulated p-value.}
#'   \item{method}{A character string indicating the type of test.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{chisq.test}}, \code{\link{lr.test}}, 
#'   \code{\link{p.test}}, \code{\link{score.test}}
#' @examples
#' # not run because the simulations take a while if B > 1000
#' ## global statistic with 4 types of record for p.test
#' #global.test(ZaragozaSeries, FUN = p.test)
#' ## global statistic with 4 types of record for chisq.test
#' #global.test(ZaragozaSeries, FUN = chisq.test)
#' ## global statistic with 4 types of record for score.test with restricted alternative
#' #global.test(ZaragozaSeries, FUN = score.test, alternative = "restricted")
#' ## global statistic with 4 types of record for lr.test with restricted alternative
#' #global.test(ZaragozaSeries, FUN = lr.test, alternative = "restricted")
#' ## global statistic with 2 types of 'almost' independent records for lr.test
#' #global.test(ZaragozaSeries, FUN = lr.test, record = c(1,0,0,1), alternative = "general")
#' 
#' @export global.test
#'
global.test <- function(X, 
                        FUN, 
                        record = c("FU" = 1, "FL" = 1, "BU" = 1, "BL" = 1), 
                        B = 1000, 
                        ...) {
  
  if (sum(record) == 0) { stop("'record' should have at least one 'TRUE' value") }

  DNAME <- deparse(substitute(X))
  METHOD <- paste("Test with global statistic for '", deparse(substitute(FUN)), 
                  "' with simulated p-value (based on ", B, " replicates)", sep = "")
  
  if (all(as.logical(record))) { # all records
    FUN_B <- function(x) {
      xrev <- series_rev(x)
      stat <- FUN(X = x, record = "upper", B = 0, ...)$statistic +
        FUN(X = x, record = "lower", B = 0, ...)$statistic +
        FUN(X = xrev, record = "upper", B = 0, ...)$statistic +
        FUN(X = xrev, record = "lower", B = 0, ...)$statistic
      return(stat)
    }
  } else if (record[1] && record[2] && !record[3] && !record[4]) { # forward records
    FUN_B <- function(x) {
      stat <- FUN(x, record = "upper", B = 0, ...)$statistic +
        FUN(x, record = "lower", B = 0, ...)$statistic
      return(stat)
    }
  } else if (!record[1] && !record[2] && record[3] && record[4]) { # backward records
    FUN_B <- function(x) {
      xrev <- series_rev(x)
      stat <- FUN(xrev, record = "upper", B = 0, ...)$statistic +
        FUN(xrev, record = "lower", B = 0, ...)$statistic
      return(stat)
    }
  } else if (record[1] && !record[2] && record[3] && !record[4]) { # upper records
    FUN_B <- function(x) {
      xrev <- series_rev(x)
      stat <- FUN(x, record = "upper", B = 0, ...)$statistic +
        FUN(xrev, record = "upper", B = 0, ...)$statistic
      return(stat)
    }
  } else if (!record[1] && record[2] && !record[3] && record[4]) { # lower records
    FUN_B <- function(x) {
      xrev <- series_rev(x)
      stat <- FUN(x, record = "lower", B = 0, ...)$statistic +
        FUN(xrev, record = "lower", B = 0, ...)$statistic
      return(stat)
    }
  } else { # others
    FUN_B <- function(x) {
      xrev <- series_rev(x)
      stat <- 0
      if (record[1]) { stat <- stat + FUN(x, record = "upper", B = 0, ...)$statistic }
      if (record[2]) { stat <- stat + FUN(x, record = "lower", B = 0, ...)$statistic }
      if (record[3]) { stat <- stat + FUN(xrev, record = "upper", B = 0, ...)$statistic }
      if (record[4]) { stat <- stat + FUN(xrev, record = "lower", B = 0, ...)$statistic }
      return(stat)
    }
  }
  
  stat <- FUN_B(X)
  
  pv <- .MonteCarlo(stat, alternative = "greater",
              FUN = FUN_B, B = B, Trows = NROW(X), Mcols = NCOL(X))
  
  names(stat) <- "Monte-Carlo"
  
  structure(list(statistic = stat, p.value = pv, 
                 method = METHOD, data.name = DNAME), class = "htest")
}
