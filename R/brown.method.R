#' @title Brown's Method on the Number of Records
#' @importFrom stats pchisq
#' @description Performs Brown's method on the p-values of \code{\link{N.test}}
#'   as proposed by Cebrián, Castillo-Mateo and Asín (2022). The null 
#'   hypothesis of the classical record model (i.e., of IID continuous RVs) is 
#'   tested against the alternative hypothesis.
#' @details 
#'   The test is implemented as given by Cebrián, 
#'   Castillo-Mateo and Asín (2022), where the p-values 
#'   \eqn{p^{(FU)}}, \eqn{p^{(FL)}}, \eqn{p^{(BU)}}, and \eqn{p^{(BL)}}
#'   of the test \code{\link{N.test}} for the four types of record are used for
#'   the statistic:
#'   \deqn{-2 \left(\log(p^{(FU)}) + \log(p^{(FL)}) + \log(p^{(BU)}) + \log(p^{(BL)})\right).}
#'   Any other combination of p-values for the test is also allowed (see 
#'   argument \code{record}).
#'   
#'   According to Brown's method (Brown, 1975) for the union of dependent 
#'   p-values, the statistic follows a \eqn{c \chi^2_{df}} distribution, 
#'   with a scale parameter \eqn{c} and \eqn{df} degrees of freedom that 
#'   depend on the covariance of the p-values. This covariances are 
#'   approximated according to Kost and McDermott (2002):
#'   \deqn{\textrm{COV}\left(-2 \log(p^{(i)}), -2 \log(p^{(j)})\right) \approx 3.263 \rho_{ij} + 0.710 \rho_{ij}^2 + 0.027 \rho_{ij}^3,}
#'   where \eqn{\rho_{ij}} is the correlation between their respective 
#'   statistics.
#'   
#'   Power studies indicate that this and \code{\link{foster.test}} using all
#'   four types of record and linear weights are the two most powerful records
#'   tests for trend detection against a linear drift model. In particular, 
#'   this test is more powerful than Mann-Kendall test against alternatives 
#'   with a linear drift in location in series of generalised Pareto variables
#'   and some cases of the generalised extreme value variables (see Cebrián, 
#'   Castillo-Mateo and Asín, 2022).
#'   
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param record Vector of length four. Each element is a logical indicating if
#'   the p-value of the test for forward upper, forward lower, backward upper 
#'   and backward lower are going to be used, respectively. Logical values or 
#'   0,1 values are accepted.
#' @param alternative Vector of length four. Each element is one of 
#'   \code{"greater"} or \code{"less"} indicating the alternative hypothesis 
#'   in every test (for forward upper, forward lower, backward upper and 
#'   backward lower records, respectively). Under the alternative hypothesis 
#'   of linear trend the FU and BL records will be greater and the FL and BU
#'   records will be less than under the null, but other combinations (e.g., 
#'   for trend in variation) could be considered.
#' @param correct Logical. Indicates, whether a continuity correction 
#'   should be applied in \code{\link{N.test}}; defaults to \code{TRUE}.
#' @return A \code{"htest"} object with elements:
#'   \item{statistic}{Value of the chi-square statistic (not scaled).}
#'   \item{parameter}{Degrees of freedom \eqn{df} and scale parameter \eqn{c}.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{fisher.method}}, \code{\link{foster.test}}, 
#'   \code{\link{N.test}}
#' @references
#' Brown M (1975). “A Method for Combining Non-Independent, One-Sided Tests of Significance.” 
#' \emph{Biometrics}. \strong{31}(4), 987–992. 
#' 
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2): 313-330. 
#' \doi{10.1007/s00477-021-02122-w}
#' 
#' Kost JT, McDermott MP (2002). “Combining Dependent P-Values.”
#' \emph{Statistics & Probability Letters}, \strong{60}(2), 183-190.
#' 
#' @examples
#' brown.method(ZaragozaSeries)
#' brown.method(ZaragozaSeries, weights = function(t) t-1)
#' brown.method(ZaragozaSeries, weights = function(t) t-1, correct = FALSE)
#' 
#' # Join p-values of upper records
#' brown.method(ZaragozaSeries, weights = function(t) t-1, record = c(1,0,1,0))
#' # Join p-values of lower records
#' brown.method(ZaragozaSeries, weights = function(t) t-1, record = c(0,1,0,1))
#' 
#' @export brown.method

brown.method <- function(X, 
                         weights = function(t) 1, 
                         record = c("FU" = 1, "FL" = 1, "BU" = 1, "BL" = 1),
                         alternative = c("FU" = "greater", "FL" = "less", "BU" = "less", "BL" = "greater"),
                         correct = TRUE) {
  
  DNAME <- deparse(substitute(X))
  X <- as.matrix(X)
  Trows <- nrow(X)
  Mcols <- ncol(X)
  
  t <- 1:Trows
  w <- weights(t)
  fun   <- deparse(weights)[2]
  if (all(w == 1)) {
    METHOD <- paste("Brown's method on the number of records")
  } else {
    METHOD <- paste("Brown's method on the weighted number of records with weights =", fun)
  }
  
  VAR <- .var_N(t, w)
  if (record[1] && record[2]) { p_NFU_NFL <- .cov_NFU_NFL(t, w) / VAR }
  if (record[1] && record[3]) { p_NFU_NBU <- .cov_NFU_NBU(t, w) / VAR }
  if (record[1] && record[4]) { p_NFU_NBL <- .cov_NFU_NBL(t, w) / VAR }
  if (record[2] && record[3]) { p_NFL_NBU <- ifelse(alternative[2] == alternative[3], 1, -1) * ifelse(record[1] && record[4], p_NFU_NBL, .cov_NFU_NBL(t, w) / VAR) }
  if (record[2] && record[4]) { p_NFL_NBL <- ifelse(alternative[2] == alternative[4], 1, -1) * ifelse(record[1] && record[3], p_NFU_NBU, .cov_NFU_NBU(t, w) / VAR) }
  if (record[3] && record[4]) { p_NBU_NBL <- ifelse(alternative[3] == alternative[4], 1, -1) * ifelse(record[1] && record[2], p_NFU_NFL, .cov_NFU_NFL(t, w) / VAR) }
  
  Cov <- function(p) {
    return(3.263 * p + 0.710 * p^2 + 0.027 * p^3)
  }
  
  k <- sum(as.logical(record))
  if (k == 0) { stop("'record' should have at least one 'TRUE' value") }
  
  E <- 2 * k
  Var <- 4 * k
  if (record[1] && record[2]) { Var <- Var + 2 * Cov(ifelse(alternative[1] == alternative[2], 1, -1) * p_NFU_NFL) }
  if (record[1] && record[3]) { Var <- Var + 2 * Cov(ifelse(alternative[1] == alternative[3], 1, -1) * p_NFU_NBU) }
  if (record[1] && record[4]) { Var <- Var + 2 * Cov(ifelse(alternative[1] == alternative[4], 1, -1) * p_NFU_NBL) }
  if (record[2] && record[3]) { Var <- Var + 2 * Cov(p_NFL_NBU) }
  if (record[2] && record[4]) { Var <- Var + 2 * Cov(p_NFL_NBL) }
  if (record[3] && record[4]) { Var <- Var + 2 * Cov(p_NBU_NBL) }
  
  cc <- Var / (2 * E)
  ff <- 2 * E^2 / Var
  
  Xrev <- series_rev(X)
  X0 <- 0
  if (record[1]) { X0 <- X0 + log(N.test(X, record = "upper", alternative = alternative[1], weights = weights, correct = correct)$p.value) }
  if (record[2]) { X0 <- X0 + log(N.test(X, record = "lower", alternative = alternative[2], weights = weights, correct = correct)$p.value) }
  if (record[3]) { X0 <- X0 + log(N.test(Xrev, record = "upper", alternative = alternative[3], weights = weights, correct = correct)$p.value) }
  if (record[4]) { X0 <- X0 + log(N.test(Xrev, record = "lower", alternative = alternative[4], weights = weights, correct = correct)$p.value) }
  X0 <- -2 * X0 / cc
  
  pv <- stats::pchisq(q = X0, df = ff, lower.tail = FALSE)
  
  names(X0) <- "X-squared"
  names(cc) <- "c"
  names(ff) <- "df"
  
  structure(list(statistic = X0, parameter = c(ff, cc),
                 p.value = pv, 
                 method = METHOD, data.name = DNAME), class = "htest")
}

.var_N <- function(v.t, v.w) {
  if (length(v.w) == 1) { VAR <- v.w^2 * sum(    (v.t[-1] - 1) / v.t[-1]^2) }
  else                  { VAR <- sum(v.w[-1]^2 * (v.t[-1] - 1) / v.t[-1]^2) }
  return(VAR)
}

.cov_NFU_NFL <- function(v.t, v.w) { 
  if (length(v.w) == 1) { COV <- - v.w^2 * sum(1 / v.t[-1]^2) }
  else                  { COV <- - sum(v.w[-1]^2 / v.t[-1]^2) }
  return(COV)
}

.cov_NFU_NBU <- function(v.t, v.w) {
  Trows <- length(v.t)
  if (length(v.w) == 1) { COV <- - v.w^2 * sum(1 / v.t[-1]^2) } # if w = 1, .cov_NFU_NBU == .cov_NFU_NFL
  else {
    x <- cumsum(rev(v.w / v.t))
    COV <- sum(v.w[-1] * rev(v.w[-Trows])) / Trows - sum(v.w[-1] / v.t[-1] * x[-1])
  }
  return(COV)
}

.cov_NFU_NBL <- function(v.t, v.w) {
  Trows <- length(v.t)
  if (length(v.w) == 1) { v.w <- rep(v.w, Trows) }
  x1 <- v.w / choose(Trows, v.t)
  x2 <- cumsum(rev(v.w / v.t))
  x3 <- 0
  lfactT <- lfactorial(Trows)
  aux <- outer(2:Trows, Trows:2, FUN = function(x, y) v.w[y] / x / (x-Trows+y-1))
  for (t in v.t[-1]) {
    x3 <- x3 + v.w[t] * sum(sweep(aux, MARGIN = 1, STATS = exp(lfactorial(v.t[-1]) + lfactorial(Trows - t) - lfactorial(v.t[-1] - t) - lfactT), FUN = `*`)[t:Trows - 1, 2:t - 1])
  }
  COV <- sum(x1[-1] * rev(v.w[-Trows]) / v.t[-1]) - sum(v.w[-1] / v.t[-1] * x2[-1]) + x3
  return(COV)
}
