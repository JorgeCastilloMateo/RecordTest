#' @title LM and LR tests with Monte Carlo global statistic
#' @importFrom stats runif
#' @description This function performs a more powerful generalization 
#'   of the Lagrange multiplier and likelihood ratio tests on records
#'   by means of the sum of the statistics of upper and lower records 
#'   and forward and backward series to study the hypothesis of the 
#'   classical record model.
#' @details Here, the Monte Carlo statistics with general alternative 
#'   \eqn{\mathcal{LM}} in \code{\link{L_lm.test}} or \eqn{\mathcal{LR}} in 
#'   \code{\link{L_lr.test}} if \code{test = 'LM'} or \code{test = 'LR'} 
#'   respectively, are joined for upper and lower records and forward and 
#'   backward series as
#'   \deqn{\mathcal{LM}^{G1} = \mathcal{LM} - \mathcal{LM}^{L} - \mathcal{LM}^{rev} + \mathcal{LM}^{L,rev}}
#'   or
#'   \deqn{\mathcal{LM}^{G2} = \mathcal{LM} + \mathcal{LM}^{L,rev}},
#'   where the superscripts \eqn{L} and \eqn{rev} means lower records and 
#'   reversed (or backwards) series, respectively. Equivalently for 
#'   \eqn{\mathcal{LR}}.
#'   
#'   The distribution of this generalized statistics is unknown, 
#'   but it can be estimated with a Monte Carlo approach.
#'
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param test A character string indicating the type of test to perform. 
#'   \code{"LM"} for \code{\link{L_lm.test}} or \code{"LR"} for 
#'   \code{\link{L_lr.test}}.
#' @param statistic A character string indicating the type of generalization
#'   to perform. \code{"G1"} or \code{"G2"} (see Details).
#' @param trend A character string indicating the type of alternative 
#'   hypothesis, positive trend in location \code{"positive"} or
#'   negative trend in location \code{"negative"}.
#' @param B An integer specifying the number of replicates used in the 
#'   Monte Carlo approach.
#' @param rdist function that simulates continuous random variables, 
#'   e.g., \code{\link{runif}} (fastest in \code{stats} package), 
#'   \code{\link{rnorm}} or \code{\link{rexp}}.
#' @param parallel If \code{TRUE}, then the Monte Carlo algorithm is done in 
#'   parallel. This can give a significant speedup on multicore machines.
#' @param numCores Allows the user to specify the amount of parallel processes 
#'   to be used if \code{parallel = TRUE}. If \code{NULL}, then the number of 
#'   logical cores is automatically detected and all available cores are used.
#' @param seed A vector of the same length as the value of \code{samples}. 
#'   Allows the user to specify the seed of each Monte Carlo simulation. 
#'   If \code{NULL}, then no seed is specified.
#' @return A list of class \code{"htest"}  with the following elements:
#'   \item{statistic}{Value of the  statistic.}
#'   \item{parameter}{Number of replicates used in the Monte Carlo test.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L_lm.test}}, \code{\link{L_lr.test}}
#' @examples
#' L_global.test(ZaragozaSeries, test = 'LM', B = 200)
#' 
#' L_global.test(ZaragozaSeries, test = 'LR', B = 200)
#' 
#' @export L_global.test

L_global.test <- function(XM_T, 
                    test = c('LM', 'LR'),
                    statistic = c('G1', 'G2'), 
                    trend = c('positive', 'negative'),
                    B = 1000, 
                    rdist = stats::runif, 
                    parallel = FALSE, 
                    numCores = 2, 
                    seed = NULL) {
  
  test <- match.arg(test)
  statistic <- match.arg(statistic)
  trend <- match.arg(trend)
  
  METHOD <- paste(test, statistic, 'record test for', trend, 'trend')
  DNAME  <- deparse(substitute(XM_T))
  XM_T   <- as.matrix(XM_T)
  Trows  <- nrow(XM_T)
  Mcols  <- ncol(XM_T)
  
  if (test == 'LM') {
    
    t <- 2:Trows
    
    LMR.fun <- function(XM_T, record) {
      
      I <- I.record(XM_T, record = record)[-1, , drop = FALSE]
      I <- sweep((sweep(I, MARGIN = 1, t, `*`) - 1)^2, MARGIN = 1, t - 1, `/`)
      
      return(sum(I))
    }
    
  } else { # test == 'LR'
    
    LMR.fun <- function(XM_T, record) { 
    
      I <- I.record(XM_T, record = record)[-1, , drop = FALSE]
    
      L0 <- rep(0, Mcols)
    
      for (i in 1:Mcols) {
      
        L0[i] <- 1 / prod(c(which(I[,i] == 1), Trows))
      }
    
      return(-2 * sum(log(L0)))
    }
  }
  
  if (statistic == 'G1') {
    
    M.fun <- function(XM_T) {
      
      XM_Trev <- apply(XM_T, 2, rev)
      
      I     <- LMR.fun(XM_T,    record = 'upper')
      IL    <- LMR.fun(XM_T,    record = 'lower')
      Irev  <- LMR.fun(XM_Trev, record = 'upper')
      ILrev <- LMR.fun(XM_Trev, record = 'lower')
      
      return(sum(I - IL - Irev + ILrev))
    }
    
  } else { # statistic == 'G2'
    
    M.fun <- function(XM_T) {
      
      XM_Trev <- apply(XM_T, 2, rev)
      
      I     <- LMR.fun(XM_T,    record = 'upper')
      ILrev <- LMR.fun(XM_Trev, record = 'lower')
      
      return(sum(I + ILrev))
    }
  }
  
  M0 <- M.fun(XM_T)
  
  #if (M0 == Inf) {
    
  #  LR.fun <- function(XM_T, record){ 
      
  #    I <- I.record(XM_T, record = record)[-1, , drop = FALSE]
      
  #    L0 <- 0
      
  #    for (i in 1:Mcols){
        
  #      L0 <- L0 + sum(log(1 / c(which(I[,i] == 1), Trows)))
  #    }
      
  #    return( -2 * L0 )
  #  }
    
  #  M0 <- M.fun(XM_T)
  #}
  
  pvalue <- MonteCarlo(M0, FUN = M.fun, trend = trend, B = B, 
                       rdist = rdist, Trows = Trows, Mcols = Mcols,
                       parallel = parallel, numCores = numCores, seed = seed)
  
  names(M0) <- "Monte-Carlo"
  names(B) <- "B"
  
  structure(list(statistic = M0, parameter = B, p.value = pvalue,
                 method = METHOD, data.name = DNAME), class = 'htest')
  
}
