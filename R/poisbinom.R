#' @name Poisson-Binomial
#' @title The Poisson Binomial Distribution
#' @aliases dpoisbinom qpoisbinom ppoisbinom rpoisbinom
#' @importFrom stats fft dbinom rbinom
#' @description Density, distribution function, quantile function and random
#'   generation for the Poisson binomial distribution with parameters 
#'   \code{size} and \code{prob}.
#'   
#'   This is conventionally interpreted as the number of successes in 
#'   \code{size * length(prob)} trials with success probabilities \code{prob}.
#' @details The Poisson binomial distribution with \code{size = 1} and 
#'   \code{prob} \eqn{= (p_1,p_2,\ldots,p_n)} has density
#'   \deqn{p(x) = \sum_{A \in F_x} \prod_{i \in A} p_i \prod_{j \in A^c} (1-p_j)}
#'   for \eqn{x=0,1,\ldots,n}; where \eqn{F_x} is the set of all subsets of 
#'   \eqn{x} integers that can be selected from \eqn{\{1,2,\ldots,n\}}.
#'   
#'   \eqn{p(x)} is computed using Hong (2013) algorithm, see the reference 
#'   below. 
#'   
#'   The quantile is defined as the smallest value \eqn{x} such that 
#'   \eqn{F(x) \ge p}, where \eqn{F} is the cumulative distribution function.
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param size The Poisson binomial distribution has \code{size} times the
#'   vector of probabilities \code{prob}.
#' @param prob Vector with the probabilities of success on each trial.
#' @param log,log.p Logical. If \code{TRUE}, probabilities \eqn{p} are given as 
#'   \eqn{\log(p)}.
#' @param lower.tail Logical. If \code{TRUE} (default), probabilities are
#'   \eqn{P(X \le x)}, otherwise, \eqn{P(X > x)}. 
#' @return 
#'   \code{dpoisbinom} gives the density, \code{ppoisbinom} gives the 
#'   distribution function, \code{qpoisbinom} gives the quantile function
#'   and \code{rpoisbinom} generates random deviates.
#'   
#'   The length of the result is determined by \code{x}, \code{q}, \code{p}
#'   or \code{n}.
#' @author Jorge Castillo-Mateo 
#' @references 
#' Hong Y (2013). 
#' “On Computing the Distribution Function for the Poisson Binomial Distribution.”
#' \emph{Computational Statistics & Data Analysis}, \strong{59}(1), 41-51.
#' \href{https://doi.org/10.1016/j.csda.2012.10.006}{doi:10.1016/j.csda.2012.10.006}
#' @rdname dpoisbinom
#' @export dpoisbinom
dpoisbinom <- function(x, size = 1, prob, log = FALSE) {
  
  if (size <= 0 | size %% 1 != 0) { stop("'size' should be a positive integer")}
  if (any(prob < 0) | any(prob > 1)) { stop("'prob' should take values in [0,1]" )}
  
  size <- rep(size, length(prob))
  prob <- rep(prob, size)
  n <- sum(size)
  
  x.0 <- which(x < 0 | x > n | x %% 1 != 0)
  
  P <- .poisbinom(n = n, prob = prob)[x + 1]
  P[x.0] <- 0
  
  if(log) { P <- log(P) } 
  
  return(P)
}

#' @rdname dpoisbinom
#' @export qpoisbinom
ppoisbinom <- function(q, size = 1, prob, lower.tail = TRUE, log.p = FALSE) {
  
  if (size <= 0 | size %% 1 != 0) { stop("'size' should be a positive integer")}
  if (any(prob < 0) | any(prob > 1)) { stop("'prob' should take values in [0,1]" )}
  
  size <- rep(size, length(prob))
  prob <- rep(prob, size)
  n <- sum(size)
  
  q.0 <- which(q < 0)
  q[q.0] <- NA
  
  P <- cumsum(.poisbinom(n = n, prob = prob))[q + 1]
  P[q.0] <- 0
  P[is.na(P)] <- 1
  
  if(log.p) { P <- log(P) }
  
  if (!lower.tail) { P <- 1 - P }
  
  return(P)
}

#' @rdname dpoisbinom
#' @export qpoisbinom
qpoisbinom <- function(p, size = 1, prob, lower.tail = TRUE, log.p = FALSE) {
  
  if (log.p) { p <- exp(p) }
  
  if (size <= 0 | size %% 1 != 0) { stop("'size' should be a positive integer")}
  if (any(prob < 0) | any(prob > 1)) { stop("'prob' should take values in [0,1]" )}
  
  size <- rep(size, length(prob))
  prob <- rep(prob, size)
  n <- sum(size)

  if (!lower.tail) { p <- 1 - p }
  
  Q <- rep(NA, length(p))
  cdf <- cumsum(.poisbinom(n = n, prob = prob))
  for (i in seq_along(p)) {
    if (p[i] < 0 | p[i] > 1) { 
      Q[i] <- NaN
    } else {
      Q[i] <- min(which(p[i] <= cdf)) - 1
    }
  }
  
  return(Q)
}

#' @rdname dpoisbinom
#' @export rpoisbinom
rpoisbinom <- function(n, size = 1, prob) {
  
  if (size <= 0 | size %% 1 != 0) { stop("'size' should be a positive integer")}
  if (any(prob < 0) | any(prob > 1)) { stop("'prob' should take values in [0,1]" )}
  
  N <- length(prob)
  
  draw <- rep(NA, n)
  for (i in 1:n) {
    draw[i] <- sum(stats::rbinom(N, size = size, prob = prob))
  }
  
  return(draw)
}

.poisbinom <- function(n, prob) {
  
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
  
  P <- abs(Re(stats::fft(c(1, a + 1i * b) / (n + 1))))
  
  return(P)
}
