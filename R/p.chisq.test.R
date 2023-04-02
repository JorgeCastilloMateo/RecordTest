#' @title Pearson's Chi-Square Test for Probabilities of Record
#' 
#' @importFrom stats pchisq
#' 
#' @description This function performs a chi-square goodness-of-fit test
#'   based on the record probabiliteis \eqn{p_t} to study the hypothesis
#'   of the classical record model (i.e., of IID continuous RVs).
#'   
#' @details 
#'   The null hypothesis of this chi-square test is that in every vector 
#'   (columns of the matrix \code{X}), the probability of record at 
#'   time \eqn{t} is \eqn{1/t} as in the classical record model, 
#'   and the alternative that the probabilities are not equal to those values. 
#'   First, the chi-square goodness-of-fit statistics to study the  null 
#'   hypothesis \eqn{H_0:\,p_t = 1/t} are calculated for each time 
#'   \eqn{t=2,\ldots,T}, where the observed value is the number of records at 
#'   time \eqn{t} in the \eqn{M} vectors and the expected value under the null
#'   is \eqn{M / t}. The test statistic is the sum of the previous \eqn{T-1} 
#'   statistics and its distribution under the null 
#'   is approximately \eqn{\chi^2_{T-1}}.
#'
#'   The chi-square approximation may not be valid with low \eqn{M}, since it
#'   requires expected values \eqn{> 5} or up to \eqn{20\%} of the expected 
#'   values are between 1 and 5. If this condition is not satisfied, a warning 
#'   is displayed. In order to avoid this problem, a \code{simulate.p.value}
#'   can be made by means of Monte Carlo simulations.
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be 
#'   calculated, "upper" or "lower".
#' @param simulate.p.value Logical. Indicates whether to compute p-values by
#'   Monte Carlo simulation. It is recommended if the function returns a 
#'   warning (see Details).    
#' @param B If \code{simulate.p.value = TRUE}, an integer specifying the 
#'   number of replicates used in the Monte Carlo estimation.
#' @return A \code{"htest"} object  with elements:
#'   \item{statistic}{Value of the chi-squared statistic.}
#'   \item{df}{Degrees of freedom.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{global.test}}, \code{\link{score.test}},
#'   \code{\link{p.record}}, \code{\link{p.regression.test}}, 
#'   \code{\link{lr.test}}
#' @references 
#' Benestad RE (2003). 
#' “How Often Can We Expect a Record Event?” 
#' \emph{Climate Research}, \strong{25}(1), 3-13.
#' \doi{10.3354/cr025003}.
#' 
#' Benestad RE (2004). 
#' “Record-Values, Nonstationarity Tests and Extreme Value Distributions.” 
#' \emph{Global and Planetary Change}, \strong{44}(1-4), 11–26. 
#' \doi{10.1016/j.gloplacha.2004.06.002}.
#' 
#' @export p.chisq.test
#' @examples
#' # Warning, M = 76 small for the value of T = 70
#' p.chisq.test(ZaragozaSeries)
#' # Simulate p-value
#' p.chisq.test(ZaragozaSeries, simulate.p.value = TRUE, B = 10000)
#'
p.chisq.test <- function(X, 
                         record = c("upper", "lower"), 
                         simulate.p.value = FALSE, 
                         B = 1000) {
  
  record <- match.arg(record)
  DNAME <- deparse(substitute(X))
  METHOD <- paste("Chi-square test on the", record, "records probabilities")

  Trows <- NROW(X)
  Mcols <- NCOL(X)
  t <- 2:Trows
  Trows_ <- Trows - 1
  
  E.R <- Mcols / t
  den <- E.R * (1 - 1 / t)
  chi.fun <- function(S) { return(sum((S - E.R)^2 / den)) }
  
  chi <- chi.fun(rowSums(.I.record(X, record = record, Trows = Trows))[-1])
  
  if (simulate.p.value) {
    METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
    SB <- matrix(stats::rbinom(n = Trows_ * B, size = Mcols, prob = 1 / t), ncol = B)
    chiB <- apply(SB, 2, chi.fun)
    pv <- sum(chiB >= chi) / B
  } else {
    pv <- stats::pchisq(chi, df = Trows_, lower.tail = FALSE)
  }

  if (((Mcols / Trows < 5 && sum(Mcols / 2:Trows < 5) > 0.2 * Trows_) ||
      Mcols / Trows < 1 || Mcols < 30) && !simulate.p.value )
    warning("Chi-squared approximation may be incorrect")

  names(chi) <- "X-squared"
  names(Trows_) <- "df"

  structure(list(statistic = chi, parameter = Trows_,
                 p.value = pv, method = METHOD, 
                 data.name = DNAME), class = "htest")
}
