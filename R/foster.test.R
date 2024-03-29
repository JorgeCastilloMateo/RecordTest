#' @title Foster-Stuart and Diersen-Trenkler Tests
#' @importFrom stats pnorm pt sd
#' @description Performs Foster-Stuart, Diersen-Trenkler and 
#'   Cebrián-Castillo-Asín records tests for trend in location, variation or 
#'   the tails. The hypothesis of the classical record model (i.e., 
#'   of IID continuous RVs) is tested against the alternative hypothesis.
#' @details 
#'   In this function, the tests are implemented as given by Foster and Stuart
#'   (1954), Diersen and Trenkler (1996, 2001) and some modifications in the
#'   standardisation of the previous statistics given by Cebrián, 
#'   Castillo-Mateo and Asín (2022). The null hypothesis is that the data come
#'   from a population with independent and identically distributed
#'   realisations. The one-sided alternative hypothesis is that the chosen
#'   statistic is greater (or less) than under the null hypothesis. The
#'   different statistics are calculated according to:
#'   
#'   If \code{statistic == "d"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} - I_{tm}^{(FL)}\right).}
#'
#'   If \code{statistic == "D"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} - I_{tm}^{(FL)} - I_{tm}^{(BU)} + I_{tm}^{(BL)}\right).}
#'
#'   If \code{statistic == "s"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} + I_{tm}^{(FL)}\right).}
#'
#'   If \code{statistic == "S"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} + I_{tm}^{(FL)} - I_{tm}^{(BU)} - I_{tm}^{(BL)}\right).}
#' 
#'   If \code{statistic == "U"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} - I_{tm}^{(BU)}\right).}
#' 
#'   If \code{statistic == "L"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(BL)} - I_{tm}^{(FL)}\right).}
#'
#'   If \code{statistic == "W"}, 
#'   \deqn{\sum_{m=1}^{M} \sum_{t=1}^{T} \omega_t \left( I_{tm}^{(FU)} + I_{tm}^{(BL)}\right).}
#'  
#'   Where \eqn{\omega_t} are weights given to the different records
#'   according to their position in the series, \eqn{I_{tm}} are the record
#'   indicators (see \code{\link{I.record}}), and \eqn{(FU)}, \eqn{(FL)}, 
#'   \eqn{(BU)}, and \eqn{(BL)} represent forward upper, forward lower,
#'   backward upper and backward lower records, respectively. The statistics 
#'   \eqn{d}, \eqn{D} and \eqn{W} may be used for trend in location;
#'   \eqn{s} and \eqn{S} may be used for trend in variation; and \eqn{U} and
#'   \eqn{L} may be used for trend in the upper and lower tails of the 
#'   distribution respectively.
#'   
#'   The statistics, say \eqn{X}, are approximately normally distributed, with
#'   \deqn{Z = \frac{X - \mu}{\sigma},}
#'   while the mean \eqn{\mu} of the particular statistic considered is simple
#'   to calculate, its variance \eqn{\sigma^2} become a cumbersome expression
#'   and some are given by Diersen and Trenkler (2001) and all of them can be
#'   easily computed out of the expression of the covariances given by Cebrián, 
#'   Castillo-Mateo and Asín (2022). 
#'   
#'   If \code{correct = TRUE}, then a continuity correction will be employed:
#'   \deqn{Z = \frac{X \pm 0.5 - \mu}{\sigma},}
#'   with ``\eqn{-}'' if the alternative is greater and ``\eqn{+}'' if the 
#'   alternative is less. Not recommended for the statistics with \eqn{\mu=0}.
#'   
#'   When \eqn{M>1}, the expression of the variance under the null hypothesis
#'   can be substituted by the sample variance in the \eqn{M} series, 
#'   \eqn{\hat{\sigma}^2}. In this case, the statistics are asymptotically
#'   \eqn{t} distributed, which is a more robust alternative against serial 
#'   correlation.
#'   
#'   If \code{permutation.test = TRUE}, the p-value is estimated by permutation
#'   simulations. This is the only method of calculating p-values that does not
#'   require that the columns of \code{X} be independent.
#'   
#'   If \code{simulate.p.value = TRUE}, the p-value is estimated by Monte Carlo
#'   simulations. If the normal asymptotic \code{statistic} \code{"D"}, 
#'   \code{"S"} or \code{"W"} is used when the length of the 
#'   series \eqn{T} is greater than 1000 or 1500, permutations or this approach
#'   are preferable due to the computational cost of calculating the variance
#'   of the statistic under the null hypothesis. The exception is \code{"D"} 
#'   without weights, which has an alternative algorithm implemented to 
#'   calculate the variance quickly.
#'   
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t - 1} then \eqn{\omega_t = t - 1}.
#' @param statistic A character string indicating the type of statistic to be 
#'   calculated, i.e., one of \code{"D"}, \code{"d"}, \code{"S"}, \code{"s"},
#'   \code{"U"}, \code{"L"} or \code{"W"} (see Details).
#' @param distribution A character string indicating the asymptotic 
#'   distribution of the statistic, \code{"normal"} or Student's
#'   \code{"t"} distribution.
#' @param alternative A character string indicating the type of alternative 
#'   hypothesis, \code{"greater"} number of records or \code{"less"} number of
#'   records.
#' @param correct Logical. Indicates, whether a continuity correction 
#'   should be done; defaults to \code{FALSE}. No correction is done if
#'   \code{simulate.p.value = TRUE}. 
#' @param permutation.test Logical. Indicates whether to compute p-values by
#'   permutation simulation (Castillo-Mateo et al. 2023). It does not require 
#'   that the columns of \code{X} be independent. If \code{TRUE} and 
#'   \code{simulate.p.value = TRUE}, permutations take precedence and 
#'   permutations are performed.
#' @param simulate.p.value Logical. Indicates whether to compute p-values by
#'   Monte Carlo simulation. If \code{permutation.test = TRUE}, permutations
#'   take precedence and permutations are performed. 
#' @param B If \code{permutation.test = TRUE} or \code{simulate.p.value = TRUE}, 
#'   an integer specifying the number of replicates used in the permutation or
#'   Monte Carlo estimation.
#' @return A \code{"htest"} object with elements:
#'   \item{statistic}{Value of the test statistic.}
#'   \item{parameter}{(If \code{distribution = "t"}) Degrees of freedom of
#'     the \eqn{t} statistic (equal to \eqn{M-1}).}
#'   \item{p.value}{P-value.}
#'   \item{alternative}{The alternative hypothesis.}
#'   \item{estimate}{(If \code{distribution = "normal"}) A vector with the
#'     value of the statistic, \eqn{\mu} and \eqn{\sigma^2}. \eqn{\sigma^2}
#'     is \code{NA} if \code{statistic} is one of \code{"D"}, \code{"S"} or 
#'     \code{"W"} (with the exception of \code{"D"} without weights); the
#'     p-value is computed with permutations or Monte Carlo simulations; and 
#'     \eqn{T > 500}.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{foster.plot}}, \code{\link{N.plot}}, 
#'   \code{\link{N.test}}
#' @references
#' Castillo-Mateo J, Cebrián AC, Asín J (2023).
#' “Statistical Analysis of Extreme and Record-Breaking Daily Maximum Temperatures in Peninsular Spain during 1960--2021.”
#' \emph{Atmospheric Research}, \strong{293}, 106934.
#' \doi{10.1016/j.atmosres.2023.106934}.
#' 
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2), 313-330. 
#' \doi{10.1007/s00477-021-02122-w}.
#' 
#' Diersen J, Trenkler G (1996). “Records Tests for Trend in Location.”
#' \emph{Statistics}, \strong{28}(1), 1-12.
#' \doi{10.1080/02331889708802543}.
#' 
#' Diersen J, Trenkler G (2001). 
#' “Weighted Records Tests for Splitted Series of Observations.”
#' In J Kunert, G Trenkler (eds.), 
#' \emph{Mathematical Statistics with Applications in Biometry: Festschrift in Honour of Prof. Dr. Siegfried Schach}, 
#' pp. 163–178. Lohmar: Josef Eul Verlag.
#' 
#' Foster FG, Stuart A (1954). 
#' “Distribution-Free Tests in Time-Series Based on the Breaking of Records.”
#' \emph{Journal of the Royal Statistical Society B}, 
#' \strong{16}(1), 1-22.
#' \doi{10.1111/j.2517-6161.1954.tb00143.x}.
#' 
#' @examples
#' # D-statistic
#' foster.test(ZaragozaSeries)
#' # D-statistic with linear weights
#' foster.test(ZaragozaSeries, weights = function(t) t - 1)
#' # S-statistic with linear weights
#' foster.test(ZaragozaSeries, statistic = "S", weights = function(t) t - 1)
#' # D-statistic with weights and t approach
#' foster.test(ZaragozaSeries, distribution = "t", weights = function(t) t - 1)
#' # U-statistic with weights (upper tail)
#' foster.test(ZaragozaSeries, statistic = "U", weights = function(t) t - 1)
#' # L-statistic with weights (lower tail)
#' foster.test(ZaragozaSeries, statistic = "L", weights = function(t) t - 1)
#' 
#' @export foster.test
foster.test <- function(X, 
                        weights = function(t) 1,
                        statistic = c("D", "d", "S", "s", "U", "L", "W"), 
                        distribution = c("normal", "t"), 
                        alternative = c("greater", "less"), 
                        correct = FALSE,
                        permutation.test = FALSE,
                        simulate.p.value = FALSE,
                        B = 1000) {
  
  if (!is.function(weights)) stop("'weights' should be a function")
  statistic    <- match.arg(statistic)
  distribution <- match.arg(distribution)
  alternative  <- match.arg(alternative)
  
  METHOD <- switch(statistic,
                   "D" = {"Foster-Stuart D-statistic test"},
                   "d" = {"Foster-Stuart d-statistic test"},
                   "S" = {"Foster-Stuart S-statistic test"},
                   "s" = {"Foster-Stuart s-statistic test"},
                   "U" = {"Forward - backward upper records test"},
                   "L" = {"Backward - forward lower records test"},
                   "W" = {"Diersen-Trenkler W-statistic test"})
  DNAME <- deparse(substitute(X))
  
  X <- as.matrix(X)
  Trows <- nrow(X)
  Mcols <- ncol(X) 
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  if (distribution == "t" && Mcols == 1) { 
    stop("If 'distribution = \"t\"', 'NCOL(X)' should be greater than 1") 
  }
  t <- 1:Trows
  w <- weights(t)
  fun   <- deparse(weights)[2]
  if (!all(w == 1)) { METHOD <- paste(METHOD, "with weights =", fun) }
  
  # Continuity correction
  n <- ifelse(correct && !permutation.test && !simulate.p.value, 0.5, 0)
  `%+-%` <- ifelse(alternative == "greater", `-`, `+`) 
  ###########
  switch(statistic,
         "D" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows) 
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I <- I.FU - I.FL - I.BU + I.BL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- ifelse(length(w) == 1,
                                    w * sqrt(4 * Mcols * sum(((t[-1] + 1 + (t[-1] - 1) / choose(Trows, t[-1])) / t[-1]^2))),
                                    ifelse(
                                      (permutation.test || simulate.p.value) && Trows > 500,
                                      NA,
                                      sqrt(4 * Mcols * (.var_N(t, w) - .cov_NFU_NFL(t, w) - .cov_NFU_NBU(t, w) + .cov_NFU_NBL(t, w)))
                                    )
                    )
                    statS <- (stat %+-% n) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows) 
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I    <- I.FU - I.FL - I.BU + I.BL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat) %+-% (n / Mcols)) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "d" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)
                      I <- I.FU - I.FL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- sqrt(2 * Mcols * (.var_N(t, w) - .cov_NFU_NFL(t, w)))
                    statS <- (stat %+-% n) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I    <- I.FU - I.FL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat)  %+-% (n / Mcols)) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "S" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows) 
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I <- I.FU + I.FL - I.BU - I.BL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- ifelse(
                      (permutation.test || simulate.p.value) && Trows > 500,
                      NA,
                      sqrt(4 * Mcols * (.var_N(t, w) + .cov_NFU_NFL(t, w) - .cov_NFU_NBU(t, w) - .cov_NFU_NBL(t, w)))
                    )
                    statS <- (stat %+-% n) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows) 
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I    <- I.FU + I.FL - I.BU - I.BL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat)  %+-% (n / Mcols)) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "s" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)
                      I <- I.FU + I.FL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 2 * Mcols * sum(w / t)
                    sigma <- sqrt(2 * Mcols * (.var_N(t, w) + .cov_NFU_NFL(t, w)))
                    statS <- (stat %+-% n - mu) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I    <- I.FU + I.FL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 2 * sum(w / t)
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat)  %+-% (n / Mcols) - mu) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "U" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows) 
                      I <- I.FU - I.BU
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- sqrt(2 * Mcols * (.var_N(t, w) - .cov_NFU_NBU(t, w)))
                    statS <- (stat %+-% n) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.BU <- .I.record(Xrev, record = "upper", Trows = Trows)    
                      I    <- I.FU - I.BU
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat)  %+-% (n / Mcols)) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "L" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I <- I.BL - I.FL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- sqrt(2 * Mcols * (.var_N(t, w) - .cov_NFU_NBU(t, w)))
                    statS <- (stat %+-% n) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FL <- .I.record(X, record = "lower", Trows = Trows)    
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I    <- I.BL - I.FL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 0
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat) %+-% (n / Mcols)) / (sigma / sqrt(Mcols))
                  }
           )
         },
         "W" = {
           switch(distribution,
                  "normal" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I <- I.FU + I.BL
                      return(sum(w * rowSums(I)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 2 * Mcols * sum(w / t)
                    sigma <- ifelse(
                      (permutation.test || simulate.p.value) && Trows > 500,
                      NA,
                      sqrt(2 * Mcols * (.var_N(t, w) + .cov_NFU_NBL(t, w)))
                    )
                    statS <- (stat %+-% n - mu) / sigma
                  },
                  "t" = {
                    stat.fun <- function(X) {
                      Xrev <- series_rev(X)
                      I.FU <- .I.record(X, record = "upper", Trows = Trows)    
                      I.BL <- .I.record(Xrev, record = "lower", Trows = Trows) 
                      I <- I.FU + I.BL
                      return(colSums(sweep(I, MARGIN = 1, STATS = w, FUN = `*`)))
                    }
                    stat  <- stat.fun(X)
                    mu    <- 2 * sum(w / t)
                    sigma <- stats::sd(stat)
                    statS <- (mean(stat)  %+-% (n / Mcols) - mu) / (sigma / sqrt(Mcols))
                  }
           )
         }
  )
  
  switch(distribution,
         "normal" = {
           if (permutation.test) {
             METHOD <- paste(METHOD, "with permuted p-value (based on", B, "permutations)")
             pv <- .permutation(stat, alternative = alternative,
                                FUN = stat.fun, B = B, X = X,
                                Trows = Trows, Mcols = Mcols)
           } else if (simulate.p.value) {
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             pv <- .MonteCarlo(stat, alternative = alternative,
                               FUN = stat.fun, B = B,
                               Trows = Trows, Mcols = Mcols)
           } else {
             pv <- switch(alternative,
                          "greater" = {stats::pnorm(statS, lower.tail = FALSE)},
                          "less"    = {stats::pnorm(statS, lower.tail = TRUE)}
             )
           }
           
           ALTERNATIVE <- paste("true 'statistic' is", alternative, "than", format(mu))
           structure(list(statistic = c("Z" = statS), p.value = pv, 
                          alternative = ALTERNATIVE,
                          estimate = c(statistic = stat, "E" = mu, "VAR" = sigma^2),
                          method = METHOD, data.name = DNAME), class = "htest")
         },
         "t" = {
           if (permutation.test) {
             METHOD <- paste(METHOD, "with permuted p-value (based on", B, "permutations)")
             statB.fun <- function(X) {
               statB <- stat.fun(X)
               return((mean(statB) - mu) / stats::sd(statB))
             }
             pv <- .permutation(statB.fun(X), alternative = alternative,
                                FUN = statB.fun, B = B, X = X,
                                Trows = Trows, Mcols = Mcols)
           } else if (simulate.p.value) {
             METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
             statB.fun <- function(X) {
               statB <- stat.fun(X)
               return((mean(statB) - mu) / stats::sd(statB))
             }
             pv <- .MonteCarlo(statB.fun(X), alternative = alternative,
                               FUN = statB.fun, B = B,
                               Trows = Trows, Mcols = Mcols)
           } else {
             pv <- switch(alternative,
                          "greater" = {stats::pt(statS, df = Mcols - 1, lower.tail = FALSE)},
                          "less"    = {stats::pt(statS, df = Mcols - 1, lower.tail = TRUE)}
             )
           }
           
           names(Mcols) <- "df"
           ALTERNATIVE <- paste("true 't' is", alternative, "than 0")
           structure(list(statistic = c("t" = statS), parameter = Mcols - 1, p.value = pv,
                          alternative = ALTERNATIVE,
                          method = METHOD, data.name = DNAME), class = "htest")
         }
  )
}
