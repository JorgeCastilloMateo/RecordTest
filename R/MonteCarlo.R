#' @title Monte Carlo simulations
#' @importFrom parallel detectCores makeCluster parLapply stopCluster 
#'   clusterEvalQ
#' @description This function performs Monte Carlo simulations when needed,
#'   it is not developed for the use of the user.
#' @keywords internal
#' @param statistic Observed value of the statistic.
#' @param trend A character string indicating if \code{statistic} is greater
#'   under the null \code{"negative"} or not \code{"positive"}.
#' @param FUN A function that computes the statistic.
#' @param B An integer specifying the number of replicates used in the Monte
#'   Carlo approach.
#' @param rdist function that simulates continuous random variables, 
#'   e.g., \code{\link{runif}} (fastest in \code{stats} package), 
#'   \code{\link{rnorm}} or \code{\link{rexp}}.
#' @param parallel If \code{TRUE}, then the Monte Carlo algorithm is done in 
#'   parallel. This can give a significant speedup on multicore machines but
#'   only if \code{B} is approximately bigger than 5000 or \code{XM_T} is big.
#' @param numCores Allows the user to specify the amount of parallel processes 
#'   to be used if \code{parallel = TRUE}. If \code{NULL}, then the number of
#'   logical cores is automatically detected and all available cores are used.
#' @param seed A vector of the same length as the value of \code{samples}. 
#'   Allows the user to specify the seed of each Monte Carlo simulation. 
#'   If \code{NULL}, then no seed is specified.
#' @param Trows,Mcols Number of rows and columns respectively, of the desired 
#'   sample \code{rdist}
#' @return Estimated p-value.
#' @author Jorge Castillo-Mateo

MonteCarlo <- function(statistic, trend = 'positive', 
                       FUN, B = 1000, rdist, parallel = FALSE, 
                       numCores = 2, seed = NULL, Trows, Mcols) {
  
  FUN_B <- function(iter) {
    
    set.seed(seed[iter])
    
    XM_TB <- matrix(rdist(Mcols * Trows), nrow = Trows, ncol = Mcols)
    
    return(FUN(XM_TB))
  }
  
  if (parallel == FALSE) {
    
    statistic_B <- sapply(1:B, FUN_B, simplify = TRUE)
    
  } else {
    
    if (is.null(numCores)) numCores <- parallel::detectCores()
    
    cl <- parallel::makeCluster(numCores)
    
    parallel::clusterEvalQ(cl = cl, library("RecordTest"))
    
    statistic_B <- parallel::parSapply(cl = cl, X = 1:B, FUN = FUN_B, simplify = TRUE)
    
    parallel::stopCluster(cl = cl)
  }
  ###################################
  
  # estimated p-value
  if (trend == 'positive') I <- ifelse(statistic_B > statistic, 1, 0)
  else                     I <- ifelse(statistic_B < statistic, 1, 0)
  
  pvalue <- sum(I) / B
  
  return(pvalue)
}