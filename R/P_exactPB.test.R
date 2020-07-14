#' @title Exact Poisson binomial test on record probabilities
#' @importFrom stats fft dbinom
#' @description This function performs an exact test  based on the record probabiliteis \eqn{p_t} to study the hypothesis of the classical record model.
#' @details The null  hypothesis  of this likelihood ratio test is that  in all the vectors (columns of matrix \code{XM_T}), the probability of record at time \eqn{t} is \eqn{1/t}.
#' The test statistic is the  total  number of records at times \eqn{t=2, ..., T}, in  the \eqn{M} vectors. Under the null, this is the sum of \eqn{M(T-1)} independent Bernoulli variables,
#' with probabilities \eqn{p_2, ...,p_2, ..., p_T, ...p_T} with \eqn{p_t=1/t}, so that its distribution is a Poisson-Binomial.
#'
#' Only unilateral alternative hypotehesis \eqn{p_t > 1/t, t=2, ..., T}  or  \eqn{p_t <1/t, t=2, ..., T} are valid, since
#' otherwise the statistic is not able to detect deviations from the null hypothesis.
#'
#' \code{\link{N_exactPB.test}}  is the same  test, but applied to only one vector, instead of M. Note in this case this test considers the probability
#' at time \eqn{t=1} (by definition of \eqn{N_t}), but the p-value  is the same.
#' @aliases P_exactPB.test N_exactPB.test
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be calculated, "upper" or "lower".
#' @param method A character string that indicates the method by which the cdf of the 
#' Poisson-Binomial distribution is calculated and therefore the p-value. 
#' \code{dft} is the preferred (and default) method, which uses the discrete Fourier transform which algorithm is given in Hong (2013).
#' \code{butler} use the algorithm given by Butler and Stephens (2016).
#' @return A \code{"htest"} object with elements:
#' \item{statistic}{Value of the likelihood ratio statistic.}
#' \item{parameter}{Number of Bernoulli independent variables summed in the statistic.}
#' \item{p.value}{P-value.}
#' \item{method}{A character string indicating the type of test performed.}
#' \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{P_chisq.test}}, \code{\link{P_regression.test}}
#' @references 
#' Butler, K. and Stephens, M.A. (2016).
#' The Distribution of a Sum of Independent Binomial Random Variables.
#' \emph{Methodology and Computing in Applied Probability}, \strong{19}(2), 557-571.
#' doi:\href{https://doi.org/10.1007/s11009-016-9533-4}{10.1007/s11009-016-9533-4}
#' 
#' Hong, Y. (2013). 
#' On computing the distribution function for the Poisson binomial distribution.
#' \emph{Computational Statistics & Data Analysis}, \strong{59}, 41-51.
#' doi:\href{https://doi.org/10.1016/j.csda.2012.10.006}{10.1016/j.csda.2012.10.006}
#' @examples
#' P_exactPB.test(ZaragozaSeries)
#' N_exactPB.test(ZaragozaSeries[, 23])
#' @export P_exactPB.test

P_exactPB.test <- function(XM_T, record = c('upper', 'lower'),
                           method = c('mixed 1', 'mixed 2', 'dft', 'butler')) {
  
  record <- match.arg(record)
  method <- match.arg(method)
    
  DNAME <- deparse(substitute(XM_T))

  XM_T <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  Trows_ <- Trows - 1
  Mcols <- ncol(XM_T)

  MN0 <- sum(I.record(XM_T, record = record)[-1, ])
  size <- Mcols * Trows_

  if (record == 'lower') {
    METHOD <- "(lower) Record indicator's exact test"
    if (MN0 == 0) {pvalue <- 0; warning('All time-series have 0 non-trivial records')}
    else pvalue <- ppoisbinom(MN0 - 1, size = Mcols, prob = 1 / 2:Trows, method = method)
  
  } else {
    METHOD <- "Record indicator's exact test"
    pvalue <- 1 - ppoisbinom(MN0, size = Mcols, prob = 1 / 2:Trows, method = method)
  }

  names(MN0) <- "Poisson-Binomial"
  names(size) <- "size"

  structure(list(statistic = MN0, parameter = size,
                 p.value = pvalue, method = METHOD, 
                 data.name = DNAME), class = 'htest')
}


#' @rdname P_exactPB.test
#' @export N_exactPB.test
N_exactPB.test <- function(XM_T, record = c('upper', 'lower'),
                           method = c('mixed 1', 'mixed 2', 'dft', 'butler')) {
  
  record <- match.arg(record)
  method <- match.arg(method)
      
  DNAME <- deparse(substitute(XM_T))

  Trows <- length(XM_T)

  NT0 <- sum(I.record(XM_T, record = record))

  if (record=='lower') {
    METHOD <- "(lower) Record counting process' exact test"
    if (NT0 == 0) {pvalue <- 0; warning('All time-series have 0 non-trivial records')}
    else pvalue <- ppoisbinom(NT0 - 1, size = 1, prob = 1 / 1:Trows, method = method)
  }
  else{
    METHOD <- "Record counting process' exact test"
    pvalue <- 1 - ppoisbinom(NT0, size = 1, prob = 1 / 1:Trows, method = method)
  }

  names(NT0) <- "Poisson-Binomial"
  names(Trows) <- "T"

  structure(list(statistic = NT0, parameter = Trows,
                 p.value = pvalue, method = METHOD, 
                 data.name = DNAME), class = 'htest')
}


ppoisbinom <- function(k, size, prob, method) {
 
  if (method == 'mixed 1') {
    
    if (size < 4) return(ppoisbinom(k = k, size = size, prob = prob, method = 'mixed 2'))
    
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
    
    index <- iter(size); index[1] <- 1
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
    
    return(sum(S[, r2]))
  
  } else if (method == 'mixed 2') {
    
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
      
    return(sum(S[, size]))
  
  } else if (method == 'dft') {
    
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
    
    return(sum(P[1:(k + 1)]))
    
  } else {
    
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
    
    return(sum(S[, m]))
  }
}

iter <- function(n) {
  
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