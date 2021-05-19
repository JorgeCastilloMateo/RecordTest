#' @title Change-point Detection Tests Based on Records
#' @importFrom stats runif
#' @description Performs change-point detection tests based on the record
#'   occurrence. The hypothesis of the classical record model (i.e., of 
#'   randomness) is tested against the alternative hypothesis that after a 
#'   certain time the series stops being stationary.
#' @details 
#'   The test is implemented as given by Castillo-Mateo (2021). The null 
#'   hypothesis is that 
#'   \deqn{H_0: p_t = 1/t, t=1,\ldots,T,}
#'   where \eqn{p_t} is the probability of (upper and/or lower) record at time 
#'   \eqn{t} and \eqn{T} is the length of the series, which corresponds to the 
#'   setup where the data come from a population with independent and 
#'   identically distributed realizations. The two-sided alternative 
#'   hypothesis is that 
#'   \deqn{H_1: p_t = 1/t, t=1,\ldots,t_0, \qquad p_t \neq 1/t, t=t_0+1,\ldots,T,}
#'   for a certain change-point \eqn{t_0}.
#'   
#'   The variables used for the statistic are
#'   \deqn{
#'   \xi_{Tt} = \omega_t \frac{\left( \frac{1}{M}\sum_{m=1}^{M}  I_{tm} -  E(I_{tm})\right)}{\sigma_T},
#'   }
#'   where \eqn{\sigma_t^2 = \sum_{k=1}^{t} \omega_k^2 Var(I_k) / M} and
#'   \eqn{\omega_t} are weights given to the different records according
#'   to their position in the \eqn{m}-th series of length \eqn{T}. If 
#'   \code{record = "upper"} or \code{= "lower"},
#'   the \eqn{I_{tm}}'s are the upper or lower record indicators, respectively
#'   (see \code{\link{I.record}}). If \code{record = "d"} or \code{= "s"},
#'   \eqn{I_{tm}} is substituted in the expressions above by 
#'   \eqn{I_{tm}^{(FU)} - I_{tm}^{(FL)}} or 
#'   \eqn{I_{tm}^{(FU)} + I_{tm}^{(FL)}}, respectively. Let 
#'   \eqn{S_{Tt} = \sum_{k=1}^{t} \xi_{Tk}}, 
#'   \eqn{s_{Tt} = \sum_{k=1}^{t} Var(\xi_{Tk})} and define the statistic
#'   \deqn{
#'   W_T(s) = S_{Tt} + \xi_{T,t+1} \frac{s - s_{Tt}}{s_{T,t+1} - s_{Tt}},
#'   }
#'   where \eqn{s \in [s_{Tt}, s_{T,t+1}]}.
#'   
#'   The test statistic is the maximum of the absolute value of the vector
#'   \deqn{
#'   B^* = \max_{1\le t \le T} |B_T(s_{Tt})|,
#'   }
#'   where \eqn{B_T(s) = W_T(s) - s W_T(1)}, \eqn{s \in [0,1]}.
#'   
#'   The probable change-point \eqn{t_0} is located where \eqn{B^*} has its 
#'   maximum. The p-value is calculated by means of the asymptotic Kolmogorov
#'   distribution. When \eqn{\omega_t \neq 1}, the asymptotic result is not 
#'   fulfilled. In that case, the p-value can be simulated using Monte Carlo 
#'   techniques with the option \code{simulate.p.value = TRUE}.
#'   
#'   As this is an asymptotic result, it has been seen that the size and power
#'   may be a little below than expected, to correct this, any of the 
#'   continuity corrections can be used:
#'   
#'   If \code{correct = "fisher"},
#'   \deqn{B^* = - \sqrt{T} \log\left(1 - \frac{B^*}{\sqrt{T}}\right).}
#'   
#'   If \code{correct = "vrbik"},
#'   \deqn{B^* = B^* + \frac{1}{6\sqrt{T}} + \frac{B^* - 1}{4T}.}
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series. Castillo-Mateo (2021)
#'   showed that the weights that get more power for this test are 
#'   \eqn{\omega_t = \sqrt{t^2 / (t-1)}}, i.e., 
#'   \code{weights = function(t) ifelse(t == 1, 0, sqrt(t^2 / (t-1)))} if
#'   \code{record = "upper"} or \code{= "lower"}.
#'   \eqn{\omega_t = \sqrt{t}}, i.e., \code{weights = function(t) sqrt(t)} if
#'   \code{record = "d"} and \eqn{\omega_t = \sqrt{t^2 / (t-2)}}, i.e., 
#'   \code{weights = function(t) ifelse(t \%in\% 1:2, 0, sqrt(t^2 / (t-2)))} if
#'   \code{record = "s"}.
#' @param record A character string that indicates the type of statistic used.
#'   If \code{record = "upper"} or \code{= "lower"}, the \eqn{I_{tm}}'s are the
#'   upper or lower record indicators, respectively. If \code{record = "d"} or 
#'   \code{= "s"}, \eqn{I_{tm}} is substituted in the expressions above by 
#'   \eqn{I_{tm}^{(FU)} - I_{tm}^{(FL)}} or 
#'   \eqn{I_{tm}^{(FU)} + I_{tm}^{(FL)}}, respectively.
#' @param correct A character string that indicates the continuity correction
#'   in the Kolmogorov distribution made to the statistic. "fisher" (Fisher and
#'   Robbins 2019), "vrbik" for the proposal by Vrbik (2020) or "none" (the 
#'   default) if no correction is made. The former shows better 
#'   size and power, but if the value of the statistic is too large it becomes 
#'   \code{NaN} and p-value \code{NA}.
#' @param simulate.p.value Logical. Indicates whether to compute p-values by
#'   Monte Carlo simulation.   
#' @param B If \code{simulate.p.value = TRUE}, an integer specifying the 
#'   number of replicates used in the Monte Carlo estimation.   
#' @return A \code{"htest"} object with elements:
#'   \item{statistic}{Value of the test statistic.}
#'   \item{p.value}{P-value.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{estimate}{The estimated change-point time.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{records}}, \code{\link{foster.test}}
#' @references 
#' Castillo-Mateo J (2021).
#' “Distribution-Free Changepoint Detection Tests Based on the Breaking of Records.”
#' Available at \href{https://arxiv.org/abs/2105.08186}{arXiv:2105.08186}
#' 
#' Fisher TJ, Robbins MW (2019). 
#' “A Cheap Trick to Improve the Power of a Conservative Hypothesis Test.”
#' \emph{The American Statistician}, \strong{73}(3), 232-242.
#' \doi{10.1080/00031305.2017.1395364}
#' 
#' Vrbik J (2020). 
#' “Deriving CDF of Kolmogorov-Smirnov Test Statistic.”
#' \emph{Applied Mathematics}, \strong{11}(3), 227-246.
#' \doi{10.4236/am.2020.113018}
#' @examples
#' change.point(ZaragozaSeries)
#' 
#' change.point(ZaragozaSeries, record = "d", 
#'   weights = function(t) sqrt(t), simulate.p.value = TRUE)
#' 
#' test.result <- change.point(rowMeans(ZaragozaSeries))
#' test.result
#' 
#' ## Not run: Load package ggplot2 to plot the changepoint
#' #library("ggplot2")
#' #records(rowMeans(ZaragozaSeries)) + 
#' #  ggplot2::geom_vline(xintercept = test.result$estimate, color = "red")
#' 
#' @export change.point

change.point <- function(X, 
                         weights = function(t) 1,
                         record = c("upper", "lower", "d", "s"),
                         correct = c("none", "fisher", "vrbik"),
                         simulate.p.value = FALSE,
                         B = 1000) {
  
  if(!is.function(weights)) stop("'weights' should be a function")
  record  <- match.arg(record)
  correct <- match.arg(correct)

  Trows <- NROW(X)
  Mcols <- NCOL(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  t     <- 1:Trows
  w     <- weights(t)
  fun   <- deparse(weights)[2]
  DNAME <- deparse(substitute(X))
  if (all(w == 1)) {
    METHOD <- paste("Records test for single changepoint detection")
  } else {
    METHOD <- paste("Records test for single changepoint detection with weights =", fun)
  }
  
  if (record == "d") {
    S      <- cumsum(w * rowSums(.I.record(X, record = "upper", Trows = Trows) - .I.record(X, record = "lower", Trows = Trows)))
    sigma2 <- c(0, Mcols * cumsum((2 * w^2 / t)[-1]))
  } else if (record == "s") {
    S      <- cumsum(w * (rowSums(.I.record(X, record = "upper", Trows = Trows) + .I.record(X, record = "lower", Trows = Trows)) - 2 * Mcols / t))
    sigma2 <- c(0, Mcols * cumsum((2 * w^2 / t * (1 - 2 / t))[-1]))
  } else { 
    S      <- cumsum(w * (rowSums(.I.record(X, record = record, Trows = Trows)) - Mcols / t))
    sigma2 <- Mcols * cumsum((w^2 / t * (1 - 1 / t)))
  }
  
  sigmaT <- sqrt(sigma2[Trows])
  W    <- S / sigmaT
  s    <- sigma2 / sigma2[Trows]
  B.   <- W - s * W[Trows]
  Babs <- abs(B.)
  Bmax <- max(Babs)
  t0   <- t[Babs == Bmax]
  
  if (simulate.p.value) {
    METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
    Bmax.fun <- function(S) {
      W    <- S / sigmaT
      B.   <- W - s * W[Trows]
      Bmax <- max(abs(B.))
      return(Bmax)
    }
    if (record == "d") {
      SB <- matrix(nrow = Trows, ncol = B)
      for (b in 1:B) SB[, b] <- cumsum(w * rowSums(matrix(.rd(Mcols, c(0, 1 / t[-1])), ncol = Mcols)))
    } else if (record == "s") {
      SB <- apply(w * (matrix(stats::rbinom(n = Trows * B, size = Mcols, prob = c(0, 2 / t[-1])), ncol = B) - c(0, 2 * Mcols / t[-1])), 2, cumsum)
    } else { 
      SB <- apply(w * (matrix(stats::rbinom(n = Trows * B, size = Mcols, prob = 1 / t), ncol = B) - Mcols / t), 2, cumsum)
    }
    BmaxB <- apply(SB, 2, Bmax.fun)
    pv <- sum(BmaxB >= Bmax) / B
  } else {
    if      (correct == "fisher") { 
      Bmax <- - sqrt(Trows) * log(1 - Bmax / sqrt(Trows)) 
    }
    else if (correct == "vrbik")  { 
      Bmax <- Bmax + 1 / (6 * sqrt(Trows)) + (Bmax - 1) / (4 * Trows) 
    }
  
    pv <- 1 - .pK(Bmax, K = 100)
  }

  structure(list(statistic = c("Kolmogorov" = Bmax), p.value = pv, 
                 alternative = "two.sided", estimate = c("probable changepoint time" = t0),
                 method = METHOD, data.name = DNAME), class = "htest")
}

.pK <- function(x, K = 100, choose = 1) {
  
  k <- 1:K
  pv <- ifelse(choose == 1,
               sqrt(2 * pi) / x * sum(exp(-(2 * k - 1)^2 * pi^2 / (8 * x^2))),
               1 - 2 * sum((-1)^(k - 1) * exp(-2 * k^2 * x^2)))
  
  return(pv)
}

.rd <- function(n, prob) {
  
  lengthp <- length(prob)
  u <- stats::runif(n = lengthp * n)
  r <- ifelse(u < prob, 1, ifelse(u > 1 - prob, -1, 0))
  return(r)
}
