#' @title Number of Records Test
#' @importFrom stats pnorm pt rbinom sd fft dbinom
#' @description Performs tests based on the (weighted) number of records, 
#'   \eqn{N^\omega}. The hypothesis of the classical record model (i.e., of IID
#'   continuous RVs) is tested against the alternative hypothesis.
#' @details 
#'   The null hypothesis is that the data come from a population with 
#'   independent and identically distributed continuous realisations. The 
#'   one-sided alternative hypothesis is that the (weighted) number of records
#'   is greater (or less) than under the null hypothesis. The 
#'   (weighted)-number-of-records statistic is calculated according to:
#'   \deqn{N^\omega = \sum_{m=1}^M \sum_{t=1}^T \omega_t I_{tm},} 
#'   where \eqn{\omega_t} are weights given to the different records
#'   according to their position in the series and \eqn{I_{tm}} are the record
#'   indicators (see \code{\link{I.record}}).
#'   
#'   The statistic \eqn{N^\omega} is exact Poisson binomial distributed
#'   when the \eqn{\omega_t}'s only take values in \eqn{\{0,1\}}. In any case,
#'   it is also approximately normally distributed, with
#'   \deqn{Z = \frac{N^\omega - \mu}{\sigma},}
#'   where its mean and variance are
#'   \deqn{\mu = M \sum_{t=1}^T \omega_t \frac{1}{t},} 
#'   \deqn{\sigma^2 = M \sum_{t=2}^T \omega_t^2 \frac{1}{t} \left(1-\frac{1}{t}\right).} 
#'   
#'   If \code{correct = TRUE}, then a continuity correction will be employed:
#'   \deqn{Z = \frac{N^\omega \pm 0.5 - \mu}{\sigma},}
#'   with ``\eqn{-}'' if the alternative is greater and ``\eqn{+}'' if the 
#'   alternative is less.
#'   
#'   When \eqn{M>1}, the expression of the variance under the null hypothesis
#'   can be substituted by the sample variance in the \eqn{M} series, 
#'   \eqn{\hat{\sigma}^2}. In this case, the statistic \eqn{N_{S}^\omega}
#'   is asymptotically \eqn{t} distributed, which is a more robust alternative
#'   against serial correlation.
#'   
#'   If \code{simulate.p.value = TRUE}, the p-value is estimated by Monte Carlo
#'   simulations.
#'   
#'   The size of the tests is adequate for any values of \eqn{T} and \eqn{M}.
#'   Some comments and a power study are given by Cebrián, Castillo-Mateo and
#'   Asín (2021).
#'   
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param record A character string indicating the type of record to be 
#'   calculated, \code{"upper"} or \code{"lower"}.
#' @param distribution A character string indicating the asymptotic 
#'   distribution of the statistic, \code{"normal"} distribution, Student's
#'   \code{"t"}-distribution or exact \code{"poisson-binomial"} distribution.
#' @param alternative A character string indicating the type of alternative 
#'   hypothesis, \code{"greater"} number of records or \code{"less"} number of
#'   records.
#' @param method (If \code{distribution = "poisson-binomial"}.) A character 
#'   string that indicates the method by which the cdf
#'   of the Poisson binomial distribution is calculated and therefore the 
#'   p-value. \code{"mixed"} is the preferred (and default) method, it is a 
#'   more efficient combination of the later algorithms. \code{"dft"} uses the 
#'   discrete Fourier transform which algorithm is given in Hong (2013). 
#'   \code{"butler"} use the algorithm given by Butler and Stephens (2016).
#' @param correct Logical. Indicates, whether a continuity correction 
#'   should be done; defaults to \code{TRUE}. No correction is done if
#'   \code{simulate.p.value = TRUE} or \code{distribution = "poisson-binomial"}. 
#' @param simulate.p.value Logical. Indicates whether to compute p-values by
#'   Monte Carlo simulation. No simulation is done if 
#'   \code{distribution = "poisson-binomial"}.    
#' @param B If \code{simulate.p.value = TRUE}, an integer specifying the 
#'   number of replicates used in the Monte Carlo estimation.
#' @return A \code{"htest"} object with elements:
#'   \item{statistic}{Value of the test statistic.}
#'   \item{parameter}{(If \code{distribution = "t"}.) Degrees of freedom of
#'     the \eqn{t} statistic (equal to \eqn{M-1}).}
#'   \item{p.value}{P-value.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{estimate}{(If \code{distribution = "normal"}) A vector with the
#'     value of \eqn{N^\omega}, \eqn{\mu} and \eqn{\sigma^2}.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{N.record}}, \code{\link{N.plot}}, 
#'   \code{\link{foster.test}}, \code{\link{foster.plot}},
#'   \code{\link{brown.method}}
#' @references 
#' Butler K, Stephens MA (2017).
#' “The Distribution of a Sum of Independent Binomial Random Variables.”
#' \emph{Methodology and Computing in Applied Probability}, \strong{19}(2), 557-571.
#' \doi{10.1007/s11009-016-9533-4}
#' 
#' Cebrián AC, Castillo-Mateo J, Asín J (2021).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}. 
#' \doi{10.1007/s00477-021-02122-w}
#' 
#' Hong Y (2013). 
#' “On Computing the Distribution Function for the Poisson Binomial Distribution.”
#' \emph{Computational Statistics & Data Analysis}, \strong{59}(1), 41-51.
#' \doi{10.1016/j.csda.2012.10.006}
#' @examples
#' # Forward Upper records
#' N.test(ZaragozaSeries)
#' # Forward Lower records
#' N.test(ZaragozaSeries, record = "lower", alternative = "less")
#' # Forward Upper records
#' N.test(series_rev(ZaragozaSeries), alternative = "less")
#' # Forward Upper records
#' N.test(series_rev(ZaragozaSeries), record = "lower")
#' 
#' # Exact test
#' N.test(ZaragozaSeries, distribution = "poisson-binom")
#' # Exact test for records in the last decade
#' N.test(ZaragozaSeries, weights = function(t) ifelse(t < 61, 0, 1), distribution = "poisson-binom")
#' # Linear weights for a more powerful test (without continuity correction)
#' N.test(ZaragozaSeries, weights = function(t) t-1, correct = FALSE)
#' 
#' @export N.test

N.test <- function(X, 
                   weights = function(t) 1, 
                   record = c("upper", "lower"), 
                   distribution = c("normal", "t", "poisson-binomial"),
                   alternative = c("greater", "less"), 
                   correct = TRUE,
                   method = c("mixed", "dft", "butler"),
                   simulate.p.value = FALSE,
                   B = 1000) {
  
  if(!is.function(weights)) stop("'weights' should be a function")
  record       <- match.arg(record)
  distribution <- match.arg(distribution)
  alternative  <- match.arg(alternative)
  method  <- match.arg(method)
  
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  if (distribution == "t" && Mcols == 1) { 
    stop("If 'distribution = \"t\"', 'NCOL(X)' should be greater than 1") 
  }
  t <- 1:Trows
  w <- weights(t)
  fun   <- deparse(weights)[2]
  DNAME <- deparse(substitute(X))
  if (all(w == 1)) {
    METHOD <- paste("Test on the number of", record, "records")
  } else {
    METHOD <- paste("Test on the weighted number of", record, "records with weights =", fun)
  }
  
  # Continuity correction
  n <- ifelse(correct && !simulate.p.value, 0.5, 0)
  `%+-%` <- ifelse(alternative == "greater", `-`, `+`) 
  ###########
  switch(distribution,
         "normal" = {
           N.fun <- function(S) sum(w * S)
           N     <- N.fun(rowSums(.I.record(X, record = record, Trows = Trows)))
           mu    <- Mcols * sum(w / t)
           sigma <- sqrt(Mcols * sum(w^2 / t * (1 - 1 / t)))
           NS    <- (N %+-% n - mu) / sigma
           
           if (simulate.p.value) {
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             SB <- matrix(stats::rbinom(n = Trows * B, size = Mcols, prob = 1 / t), ncol = B)
             NB <- apply(SB, 2, N.fun)
             pv <- switch(alternative,
                          "greater" = {sum(NB >= N) / B},
                          "less"    = {sum(NB <= N) / B}
             )
           } else {
             pv <- switch(alternative,
                          "greater" = {stats::pnorm(NS, lower.tail = FALSE)},
                          "less"    = {stats::pnorm(NS, lower.tail = TRUE)}
             )
           }
           
           ALTERNATIVE <- paste("true 'N' is", alternative, "than", format(mu))
           structure(list(statistic = c("Z" = NS), p.value = pv, 
                          alternative = ALTERNATIVE,
                          estimate = c("N" = N, "E" = mu, "VAR" = sigma^2),
                          method = METHOD, data.name = DNAME), class = "htest")
         },
         "t" = {
           # argument S is I
           N.fun <- function(S) colSums(sweep(S, MARGIN = 1, STATS = w, FUN = `*`))
           N     <- N.fun(.I.record(X, record = record, Trows = Trows))
           mu    <- sum(w / t)
           sigma <- stats::sd(N)
           NS    <- (mean(N)  %+-% (n / Mcols)  - mu) / (sigma / sqrt(Mcols))
           
           if (simulate.p.value) {
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             IB <- array(stats::rbinom(n = Trows * Mcols * B, size = 1, prob = 1 / t), dim = c(Trows, Mcols, B))
             NB <- apply(IB, 3, N.fun)
             NSB <- (apply(NB, 2, mean) - mu) / (apply(NB, 2, stats::sd) / sqrt(Mcols))
             pv <- switch(alternative,
                          "greater" = {sum(NSB >= NS) / B},
                          "less"    = {sum(NSB <= NS) / B}
             )
           } else {
             pv <- switch(alternative,
                          "greater" = {stats::pt(NS, df = Mcols - 1, lower.tail = FALSE)},
                          "less"    = {stats::pt(NS, df = Mcols - 1, lower.tail = TRUE)}
             )
           }
           
           names(Mcols) <- "df"
           ALTERNATIVE <- paste("true 't' is", alternative, "than 0")
           structure(list(statistic = c("t" = NS), parameter = Mcols - 1, p.value = pv, 
                          alternative = ALTERNATIVE,
                          method = METHOD, data.name = DNAME), class = "htest")
         },
         "poisson-binomial" = {
           if (!all(w %in% 0:1)) {stop("If 'distribution = \"poisson-binomial\"', 'weights' should take values 0's or 1's")}
           N.fun <- function(S) sum(w * S)
           N <- N.fun(rowSums(.I.record(X, record = record, Trows = Trows)))
           if (length(w) > 1) { t <- t[which(w == 1)] }
           
           pv <- switch(alternative,
             "greater" = {.ppoisbinom(N - 1, size = Mcols, prob = 1 / t, lower.tail = FALSE, method = method)},
             "less"    = {.ppoisbinom(N, size = Mcols, prob = 1 / t, lower.tail = TRUE, method = method)}
           )
           
           ALTERNATIVE <- paste("true 'Poisson-binomial' is", alternative, "than", format(Mcols * sum(w / 1:Trows)))
           structure(list(statistic = c("Poisson-binomial" = N), p.value = pv, 
                          alternative = ALTERNATIVE,# parameter = c("M" = Mcols, "prob" = 1 / t),
                          method = METHOD, data.name = DNAME), class = "htest")
         }
  )
}

.ppoisbinom <- function(k, size, prob, lower.tail, method) {
  
  if (k < 0) { return(ifelse(lower.tail, 0, 1)) }

  if (method == "mixed") {
    
    if (size < 4) {
      return(.ppoisbinom(k = k, size = size, prob = prob, lower.tail = lower.tail, method = "mixed 2"))
    }
    
    n <- length(prob)
    L <- ceiling(n / 2)
    
    a <- rep(1, n)
    b <- rep(1, n)
    
    w <- 2 * pi / (n + 1)
    
    for (l in 1:L) {
      z <- 1 + prob * (-1 + cos(w * l) + 1i * sin(w * l)) 
      d <- exp(sum(log(abs(z))))
      a[l] <- d * cos(sum(Arg(z)))
      b[l] <- d * sin(sum(Arg(z)))
      a[n + 1 - l] <- a[l]
      b[n + 1 - l] <- -b[l]
    }
    
    x <- a + 1i * b
    P <- abs(Re(stats::fft(c(1, x) / (n + 1))))
    
    ########################################################
    
    index <- .iter(size); index[1] <- 1
    r2 <- length(index) + 1
    S <- matrix(0, nrow = k + 1, ncol = r2)
    S[1:min(n + 1, k + 1), 1] <- P[1:min(n + 1, k + 1)]
    index2 <- 1
    
    for (r in 2:r2) {
      if (index[r - 1] == 1) {  
        index2 <- index2 + 1
        for (j in 0:min(index2 * n, k)) {
          for (i in max(0, j - n):j) {
            S[j + 1, r] <- S[j + 1, r] + S[i + 1, r - 1] * P[j - i + 1]
          }
        }
      } else if (index[r - 1] == 0) {
        index2 <- 2 * index2
        for (j in 0:min(index2 * n, k)) {
          for (i in 0:j) {
            S[j + 1, r] <- S[j + 1, r] + S[i + 1, r - 1] * S[j - i + 1, r - 1]
          }
        }
      }
    }
    
    return(ifelse(lower.tail, sum(S[, r2]), 1 - sum(S[, r2])))
    
  } else if (method == "mixed 2") {
    
    n <- length(prob)
    L <- ceiling(n / 2)
    
    a <- rep(1, n)
    b <- rep(1, n)
    
    w <- 2 * pi / (n + 1)
    
    for (l in 1:L) {
      z <- 1 + prob * ( -1 + cos(w * l) + 1i * sin(w * l)) 
      d <- exp(sum(log(abs(z))))
      a[l] <- d * cos(sum(Arg(z)))
      b[l] <- d * sin(sum(Arg(z)))
      a[n + 1 - l] <- a[l]
      b[n + 1 - l] <- -b[l]
    }
    
    x <- a + 1i * b
    P <- abs(Re(stats::fft(c(1, x) / (n + 1))))
    
    ########################################################
    
    S <- matrix(0, nrow = k + 1, ncol = size)
    S[1:min(n + 1, k + 1), 1] <- P[1:min(n + 1, k + 1)]
    
    if (size > 1) {
      for (r in 2:size) {
        for (j in 0:min(r * n, k)) {
          for (i in max(0, j - n):j) {
            S[j + 1, r] <- S[j + 1, r] + S[i + 1, r - 1] * P[j - i + 1]
          }
        }
      }
    }
    
    return(ifelse(lower.tail, sum(S[, size]), 1 - sum(S[, size])))
    
  } else if (method == "dft") {
    
    size <- rep(size, length(prob))
    prob <- rep(prob, size)
    
    n <- sum(size)
    L <- ceiling(n / 2)
    
    a <- rep(1, n)
    b <- rep(1, n)
    
    w <- 2 * pi / (n + 1)
    
    for (l in 1:L) {
      z <- 1 + prob * (-1 + cos(w * l) + 1i * sin(w * l)) 
      d <- exp(sum(log(abs(z))))
      a[l] <- d * cos(sum(Arg(z)))
      b[l] <- d * sin(sum(Arg(z)))
      a[n + 1 - l] <- a[l]
      b[n + 1 - l] <- -b[l]
    }
    
    x <- a + 1i * b
    P <- abs(Re(stats::fft(c(1, x) / (n + 1))))
    
    return(ifelse(lower.tail, sum(P[1:(k + 1)]), 1 - sum(P[1:(k + 1)])))
    
  } else { # method = "butler"
    
    m <- length(prob)
    S <- matrix(0, nrow = k + 1, ncol = m)
    
    for (j in 0:min(size, k)) S[j + 1, 1] <- stats::dbinom(j, size, prob[1])
    
    for (r in 2:m) {
      for (j in 0:min(r * size, k)) {
        for (i in max(0, j - size):j) {
          S[j + 1, r] <- S[j + 1, r] + S[i + 1, r - 1] * stats::dbinom(j - i, size, prob[r])
        }
      }
    }
    
    return(ifelse(lower.tail, sum(S[, m]), 1 - sum(S[, m])))
  }
}

.iter <- function(n) {
  
  index <- c()
  
  while (n != 1) {
    if (n %% 2 == 1) {
      n <- n - 1
      index <- c(1, index)
    } 
    n <- n / 2; index <- c(0, index)
  }
  
  return(index)
}
