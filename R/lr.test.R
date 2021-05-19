#' @title Likelihood-Ratio Test for the Likelihood of the Record Indicators
#' @importFrom stats pchisq
#' @description This function performs likelihood-ratio tests
#'   for the likelihood of the record indicators \eqn{I_t} to study the 
#'   hypothesis of the classical record model.
#' @details 
#'   The null hypothesis of the likelihood-ratio tests is that in every vector 
#'   (columns of the matrix \code{X}), the probability of record at 
#'   time \eqn{t} is \eqn{1 / t} as in the classical record model (i.e., 
#'   sequences of independent and identically distributed realizations), and 
#'   the alternative depends on the \code{alternative} and \code{probabilities}
#'   arguments. The probability at time \eqn{t} is any value, but equal in the
#'   \eqn{M} series if \code{probabilities = "equal"}  or different in the 
#'   \eqn{M} series if \code{probabilities = "different"}. The alternative 
#'   hypothesis is more specific in the first case than in the second one.
#'   Furthermore, the \code{"two.sided"} \code{alternative} is tested with 
#'   the usual likelihood ratio statistic, while the one-sided 
#'   \code{alternatives} use specific statistics based on likelihoods. 
#'   (See Cebrián, Castillo-Mateo and Asín (2021) for details on these tests.)
#'
#'   If \code{alternative = "two.sided" & probabilities = "equal"}, under the
#'   null, the likelihood ratio statistic has an asymptotic \eqn{\chi^2} 
#'   distribution with \eqn{T-1} degrees of freedom. It has been seen that for
#'   the approximation to be adequate \eqn{M} must be between 4 and 5 times 
#'   greater than \eqn{T}. Otherwise, a \code{simulate.p.value} is recommended.
#'   
#'   If \code{alternative = "two.sided" & probabilities = "different"}, the 
#'   asymptotic behavior is not fulfilled, but the Monte Carlo approach to 
#'   simulate the p-value is applied. This statistic is the same as \eqn{\ell} 
#'   below multiplied by a factor of 2, so the p-value is the same.
#'   
#'   If \code{alternative} is one-sided and \code{probabilities = "equal"},
#'   the statistic of the test is
#'   \deqn{-2 \sum_{t=2}^T \left\{-S_t \log\left(\frac{tS_t}{M}\right)+(M-S_t)\left( \log\left(1-\frac{1}{t}\right) - \log\left(1-\frac{S_t}{M}\right) I_{\{S_t<M\}} \right) \right\} I_{\{S_t > M/t\}}.}
#'   The p-value of this test is estimated with Monte Carlo simulations,
#'   because the computation of its exact distribution is very expensive.   
#'   
#'   If \code{alternative} is one-sided and \code{probabilities = "different"},
#'   the statistic of the test is
#'   \deqn{\ell = \sum_{t=2}^T  S_{t} \log(t-1) - M \log\left(1-\frac{1}{t}\right).}
#'   The p-value of this test is estimated with Monte Carlo simulations. 
#'   However, it is equivalent to the statistic of the weighted number of 
#'   records \code{\link{N.test}} with weights \eqn{\omega_t = \log(t-1)} 
#'   \eqn{(t=2,\ldots,T)}.
#' 
#' @inheritParams score.test   
#' 
#' @return A list of class \code{"htest"} with the following elements:
#'   \item{statistic}{Value of the statistic.}
#'   \item{parameter}{Degrees of freedom of the approximate \eqn{\chi^2} 
#'     distribution.}
#'   \item{p.value}{(Estimated) P-value.}
#'   \item{method}{A character string indicating the type of test.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   \item{alternative}{A character string indicating the alternative
#'     hypothesis.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{global.test}}, \code{\link{score.test}}
#' @references 
#' Cebrián A, Castillo-Mateo J, Asín J (2021).
#' “Record Tests to detect non stationarity in the tails with an application to climate change.”
#' Available at Research Square \doi{10.21203/rs.3.rs-214787/v1}
#' 
#' @examples
#' set.seed(23)
#' # two-sided and different probabilities of record, always simulated the p-value
#' lr.test(ZaragozaSeries, probabilities = "different")
#' # equal probabilities
#' lr.test(ZaragozaSeries, probabilities = "equal")
#' # equal probabilities with simulated p-value
#' lr.test(ZaragozaSeries, probabilities = "equal", simulate.p.value = TRUE)
#' 
#' # one-sided and different probabilities of record
#' lr.test(ZaragozaSeries, alternative = "greater", probabilities = "different")
#' # different probabilities with simulated p-value
#' lr.test(ZaragozaSeries, alternative = "greater", probabilities = "different", 
#'   simulate.p.value = TRUE)
#' # equal probabilities, always simulated the p-value
#' lr.test(ZaragozaSeries, alternative = "greater", probabilities = "equal")
#' @export lr.test
#' 
lr.test <- function(X, 
                    record = c("upper", "lower"), 
                    alternative = c("two.sided", "greater", "less"),
                    probabilities = c("different", "equal"), 
                    simulate.p.value = FALSE,
                    B = 1000) {
  
  record <- match.arg(record)
  alternative <- match.arg(alternative)
  probabilities <- match.arg(probabilities)
  METHOD <- paste("Likelihood-ratio test for", record, "record indicators")
  DNAME <- deparse(substitute(X))
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  t <- 2:Trows
  
  if (alternative == "two.sided") {
    switch(probabilities,
           "equal" = {
             LR.fun <- function(S) {
               L0 <- (1 / t)^S * (1 - 1 / t)^(Mcols - S)
               L1 <- (S / Mcols)^S * (1 - S / Mcols)^(Mcols - S)
               LR <- 2 * (sum(log(L1)) - sum(log(L0)))
               return(LR)
             }
             
             LR0 <- LR.fun(rowSums(.I.record(X, record = record, Trows = Trows))[-1])
             
             if (simulate.p.value) {
               METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
               SB <- matrix(stats::rbinom(n = (Trows - 1) * B, size = Mcols, prob = 1 / t), ncol = B)
               LRB <- apply(SB, 2, LR.fun)
               pv <- sum(LRB >= LR0) / B
             } else {
               pv <- stats::pchisq(q = LR0, df = Trows - 1, lower.tail = FALSE)
             }
             
             names(LR0) <- "X-squared"
             names(Trows) <- "df"
             
             structure(list(statistic = LR0, parameter = Trows - 1, p.value = pv, 
                            method = METHOD, data.name = DNAME, 
                            alternative = paste(alternative, "with", probabilities, "probabilities")),
                       class = "htest")
           },
           "different" = {
             # likelihood ratio function
             logt <- log(t - 1)
             LR.fun <- function(S) { sum(S * logt) }
             ###################################
             # likelihood ratio
             LR0 <- LR.fun(rowSums(.I.record(X, record = record, Trows = Trows)[-1, , drop = FALSE]))
             ###################################
             # p-value Monte-Carlo
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             SB <- matrix(stats::rbinom(n = (Trows - 1) * B, size = Mcols, prob = 1 / t), ncol = B)
             LRB <- apply(SB, 2, LR.fun)
             pv <- sum(LRB >= LR0) / B
             
             LR0 <- 2 * (LR0 - Mcols * sum(log((t - 1) / t)))
             ###################################
             
             names(LR0) <- "LR"
             
             structure(list(statistic = LR0, p.value = pv,
                            method = METHOD, data.name = DNAME, 
                            alternative = paste(alternative, "with", probabilities, "probabilities")),
                       class = "htest")
           })
  } else { # alternative %in% c("greater", "less")
    switch(probabilities,
           "equal" = {
             # likelihood ratio function
             LR.fun <- function(S) {
               sum(ifelse(S > Mcols / t, 
                          S * log(S * t / Mcols) + 
                            ifelse(S == Mcols, 
                                   0, 
                                   (Mcols - S) * log(t * (Mcols - S) / (Mcols * (t - 1)))), 
                          0)
               )
             }
             ###################################
             # LR statistic
             LR0 <- LR.fun(rowSums(.I.record(X, record = record, Trows = Trows)[-1, , drop = FALSE]))
             ###################################
             # p-value Monte-Carlo
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             SB <- matrix(stats::rbinom(n = (Trows - 1) * B, size = Mcols, prob = 1 / t), ncol = B)
             LRB <- apply(SB, 2, LR.fun)
             pv <- switch (alternative,
                           "greater" = {sum(LRB >= LR0) / B},
                           "less"    = {sum(LRB <= LR0) / B}
             )
             ###################################
             LR0 <- 2 * LR0
             
             names(LR0) <- "LR"
             
             structure(list(statistic = LR0, p.value = pv, 
                            method = METHOD, data.name = DNAME, 
                            alternative = paste(alternative, "with", probabilities, "probabilities")),
                       class = "htest")
           },
           "different" = {
             # likelihood ratio function
             logt <- log(t - 1)
             LR.fun <- function(S) { sum(S * logt) }
             ###################################
             # LR statistic
             LR0 <- LR.fun(rowSums(.I.record(X, record = record, Trows = Trows)[-1, , drop = FALSE]))
             ###################################
             # p-value Monte-Carlo
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             SB <- matrix(stats::rbinom(n = (Trows - 1) * B, size = Mcols, prob = 1 / t), ncol = B)
             LRB <- apply(SB, 2, LR.fun)
             pv <- switch (alternative,
                           "greater" = {sum(LRB >= LR0) / B},
                           "less"    = {sum(LRB <= LR0) / B}
             )
             ###################################
             LR0 <- LR0 - Mcols * sum(log((t - 1) / t))
             
             names(LR0) <- "LR"
             
             structure(list(statistic = LR0, p.value = pv, 
                            method = METHOD, data.name = DNAME, 
                            alternative = paste(alternative, "with", probabilities, "probabilities")),
                       class = "htest")
           })
  }
}
