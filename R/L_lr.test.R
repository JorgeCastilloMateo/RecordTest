#' @title Likelihood ratio tests on record times
#' @importFrom stats pchisq
#' @description This function performs an asymptotic likelihood ratio test or 
#'   a Monte Carlo approach based on the record times \eqn{L_i} to study the 
#'   hypothesis of the classical record model.
#' @details The null hypothesis of this likelihood ratio tests is that in all
#'   the series, \eqn{m=1,\ldots,M}, the probability of record at time \eqn{t}
#'   is \eqn{1/t}, and the alternative depends on the \code{alternative} 
#'   argument. The probability at time \eqn{t} is any value, but equal in the 
#'   \eqn{M} series if \code{alternative = '='} or different in the 
#'   \eqn{M} series if \code{alternative = '!='}. The alternative hypothesis
#'   is more specific in the first case than in the second one.
#'
#'   If \code{alternative = '='}, under the null, the likelihood ratio 
#'   statistic has an asymptotic \eqn{\chi^2} distribution with \eqn{T-1} 
#'   degrees of freedom. It has been seen that for the approximation to be 
#'   adequate \eqn{M} must be between 4 and 5 times greater than \eqn{T}. 
#'   
#'   If \code{alternative = '!='}, the asymptotic behaviour is not fulfilled, 
#'   but a Monte Carlo approach can be applied.
#'   
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record, 
#'   "upper" or "lower".
#' @param alternative A character indicating if the alternative hypothesis 
#'   assume all series with the same distribution \code{"="} (asymptotic 
#'   chi-squared) or different distribution \code{"!="} (Monte Carlo 
#'   simulated p-value).
#' @param trend A character string indicating the type of alternative 
#'   hypothesis, positive trend in location \code{"positive"} or
#'   negative trend in location \code{"negative"}.
#' @param B An integer specifying the number of replicates used in the 
#'   Monte Carlo approach. Only used if \code{alternative = '!='}.
#' @param rdist A function that simulates continuous random variables, 
#'   e.g., \code{\link{runif}} (fastest in \code{stats} package), 
#'   \code{\link{rnorm}} or \code{\link{rexp}}. Only used if 
#'   \code{alternative = '!='}.
#' @param parallel If \code{TRUE}, then the Monte Carlo algorithm is done in 
#'   parallel. This can give a significant speedup on multicore machines.
#'   Only used if \code{alternative = '!='}.
#' @param numCores Allows the user to specify the amount of parallel processes
#'   to be used if \code{parallel = TRUE} and \code{alternative = '!='}. 
#'   If \code{NULL}, then the number of logical cores is automatically 
#'   detected and all available cores are used.
#' @param seed A vector of the same length as the value of \code{samples}. 
#'   Allows the user to specify the seed of each Monte Carlo simulation. 
#'   If \code{NULL}, then no seed is specified. Only used if 
#'   \code{alternative = '!='}.
#' @return A list of class \code{"htest"}  with the following elements:
#'   \item{statistic}{Value of the  statistic.}
#'   \item{parameter}{Degrees of freedom of the approximate \eqn{\chi^2} 
#'     distribution, or number of replicates used in the Monte Carlo test.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   \item{alternative}{A character string indicating one of both alternative
#'     hypothesis.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L_lm.test}}, \code{\link{L_global.test}}
#' @examples
#' L_lr.test(ZaragozaSeries, B = 200)
#' 
#' L_lr.test(series_rev(ZaragozaSeries), alternative = '=', trend = 'negative')
#' @export L_lr.test

L_lr.test <- function(XM_T, record = c('upper', 'lower'), 
                      alternative = c('!=', '='), trend = c('positive', 'negative'),
                      B = 1000, rdist = stats::runif,
                      parallel = FALSE, numCores = 2, seed = NULL) {
  
  record <- match.arg(record)
  alternative <- match.arg(alternative)
  trend <- match.arg(trend)
  METHOD <- paste("Likelihood ratio test on", record, "records for", trend, "trend")
  DNAME <- deparse(substitute(XM_T))
  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  
  if      (record == 'lower' && trend == 'positive') trend <- 'negative'
  else if (record == 'lower' && trend == 'negative') trend <- 'positive'
  
  if (alternative == '=') {
  
    XM_T <- I.record(XM_T, record = record)[-1, , drop = FALSE]
    
    # null hypothesis
    a <- list()
    L0 <- rep(0, Mcols)
    
    for (m in 1:Mcols) {
      
      # record times - 1 on series m
      a[[m]] <- which(XM_T[,m] == 1)
      
      # likelihood on series m under H_0
      L0[m] <- 1 / prod(c(a[[m]], Trows))
    }
    ###################################
    
    # alternative hypothesis
    a <- unlist(a)
    df <- Trows - 1 # degrees of freedom
    alpha <- L1 <- rep(0, df)
    
    for (t in 1:df){
      
      # number of records at time t-1 in the M series
      alpha[t] <- length(which(a == t))
      
      # components of the product sequence equal to the likelihood under H_1
      L1[t] <- (alpha[t] / Mcols)^alpha[t] * (1 - alpha[t]/Mcols)^(Mcols - alpha[t])
    }
    ###################################
    
    # likelihood ratio
    ELL0 <- 2 * ( sum(log(L1)) - sum(log(L0)) )
    
    pvalue <- stats::pchisq(ELL0, df = df, lower.tail = (trend == 'negative'))
    
    names(ELL0) <- "X-squared"
    names(df) <- "df"
    
    structure(list(statistic = ELL0, parameter = df, p.value = pvalue, 
                   method = METHOD, data.name = DNAME, 
                   alternative = "restricted"), class = 'htest')
  
  } else {
    
    # likelihood ratio function
    ELL0.fun <- function(XM_T){ 
      
      I <- I.record(XM_T, record = record)[-1, , drop = FALSE]
      
      a <- list()
      L0 <- rep(0, Mcols)
      
      for (i in 1:Mcols){
        
        a[[i]] <- which(I[,i] == 1)
        
        L0[i] <- 1 / prod(c(a[[i]], Trows))
      }
      
      return(-2 * sum(log(L0)))
    }
    ###################################
    
    # likelihood ratio
    ELL0 <- ELL0.fun(XM_T)
    ###################################
    
    pvalue <- MonteCarlo(ELL0, trend = trend, FUN = ELL0.fun, B = B,
                         rdist = rdist, parallel = parallel, numCores = numCores,
                         seed = seed, Trows = Trows, Mcols = Mcols)
    ###################################
    
    names(ELL0) <- "Monte-Carlo"
    names(B) <- "B"
    
    structure(list(statistic = ELL0, parameter = B, p.value = pvalue,
                   method = METHOD, data.name = DNAME, 
                   alternative = "general"), class = 'htest')
  }
}
