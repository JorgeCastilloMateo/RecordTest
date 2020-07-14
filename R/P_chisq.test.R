#' @title Chi-Square test on record probabilities
#' @importFrom stats pchisq
#' @description This function performs a chi-square test based on the record 
#'   probabiliteis \eqn{p_t} to study the hypothesis of the classical record 
#'   model.
#' @details The null hypothesis of this chi-square test is that in all the 
#'   vectors (columns of matrix \code{XM_T}), the probability of record at 
#'   time \eqn{t} is \eqn{1/t}, and the alternative that the probabilities 
#'   are not equal to those values. First, the chi-square goodness of fit
#'   statistics to study the  null hypotehsis \eqn{H_0:\,p_t = 1/t} are 
#'   calculated for each time \eqn{t=2,\ldots,T}, where the observed value 
#'   is the number of records  at time \eqn{t} in the \eqn{M} vectors and the
#'   expected value under the null is \eqn{M/t}. The test statistic is the sum
#'   of the previous \eqn{T-1} statistics and its distribution under the null 
#'   is approximately \eqn{\chi^2_{T-1}}.
#'
#'   The chi-square approximation may not be valid with low \eqn{M}, since it
#'   requires expected values \eqn{> 5} or up to 20 \% of the expected values
#'   is between 1 and 5. If this condition is not satisfied, a warning is 
#'   displayed. In order to avoid this problem, a \code{correction} can be 
#'   made.
#'
#' @param XM_T A matrix.
#' @param record A character string indicating the type of record to be 
#'   calculated,  "upper" or "lower".
#' @param correction Logical flag. If \code{TRUE}, a generalization of the 
#'   record indicator random variables is calculated by making all expected 
#'   values greater than 5.
#' @return A  \code{"htest"} object  with elements:
#'   \item{statistic}{Value of the chi-squared statistic.}
#'   \item{df}{Degrees of freedom of the approximate chi-squared.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{P_exactPB.test}}, \code{\link{P_regression.test}}
#' @export P_chisq.test
#' @examples
#' P_chisq.test(ZaragozaSeries)
#' P_chisq.test(ZaragozaSeries, correction = FALSE)

P_chisq.test <- function(XM_T, record = c('upper', 'lower'), correction = TRUE){
  
  METHOD <- "Chi-Square test on record probabilities"
  DNAME <- deparse(substitute(XM_T))

  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  
  if (correction) {
    
    l <- l.fun(Trows, Mcols)
    Mt <- Ml.rec(XM_T, record = record, l = l[,1])
    ER <- Mcols * l[,2]
    Trows_ <- nrow(l)
  
  } else {
    
    Mt <- M.record(XM_T, record = record)
    ER <- Mcols / 2:Trows
    Trows_ <- Trows - 1
  }

  chi <- rep(0, Trows_)
  
  ENR <- Mcols - ER
  chi <- ((Mt[-1] - ER)^2 / ER + (Mcols - Mt[-1] - ENR)^2 / ENR)
  
  chi <- sum(chi)
  pvalue <- stats::pchisq(chi, df = Trows_, lower.tail = FALSE)


  if (((Mcols / Trows < 5 && length(which(Mcols / 2:Trows < 5)) > 0.2 * Trows_) ||
      Mcols / Trows < 1 || Mcols < 30) & !correction )
    warning("Chi-squared approximation may not be valid")

  names(chi) <- "X-squared"
  names(Trows_) <- "df"

  structure(list(statistic = chi, parameter = Trows_,
                 p.value = pvalue, method = METHOD, 
                 data.name = DNAME), class = 'htest')
}

 
l.fun <- function(last, M) {
  
  if (M <= 10) stop("the number of series M is not enough")
  
  i <- 0
  l <- c()
  pt <- c()
  
  while (1 / last <= 5 / M) {
    
    i <- i + 1
    l[i] <- 0
    pt[i] <- 1 / last
    prod <- (1 - 1 / last)
    
    while (pt[i] <= 5 / M) {

      l[i] <- l[i] + 1
      prod <- prod * (1 - 1 / (last - l[i]))
      pt[i] <- 1 - prod
    }
    
    last <- last - (l[i] + 1)
  }
  
  if (last == 1) return(cbind(rev(l), rev(pt)))
  
  I <- cbind(rep(0, last - 1), 1 / 2:last)
  
  return(rbind(I, cbind(rev(l), rev(pt))))
}

Ml.rec <- function(XM_T, record, l) {
  
  XM_T <- I.record(XM_T, record = record)
  
  X <- matrix(nrow = length(l) + 1, ncol = ncol(XM_T))
  
  X[1, ] <- rep(1, ncol(XM_T))
  
  for (k in 2:(length(l) + 1)) {
    
    j <- k + sum(l[1:(k-1)])
    i <- j - l[k-1]

    if (i == j) X[k, ] <- XM_T[i, ]
    else        X[k, ] <- apply(as.matrix(XM_T[i:j, ]), 2, max)
  }
  
  X <- apply(X, 1, sum)
  
  return(X)
}