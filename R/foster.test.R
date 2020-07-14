#' @title Foster-Stuart, Diersen-Trenkler or t Tests for Randomness
#' @importFrom stats pnorm pt sd
#' @description Performs Foster-Stuart, Diersen-Trenkler 
#'   or t record tests for trend in location based on their normal 
#'   or t asymptotic distribution, respectively.
#'   The null hypothesis of randomness is tested against the 
#'   alternative hypothesis.
#' @details 
#'   In this function, the tests are implemented as given by Foster and Stuart
#'   (1954), Diersen and Trenkler (1996, 2001) and some modifications
#'   standardizing the previous statistics. Let \eqn{(I_{t,m})} be the record 
#'   indicator random variables of the mth series. The number of records up to
#'   time \eqn{T} of the mth series, \eqn{N_{T,m} = \sum_{t=2}^{T} I_{t,m}}, 
#'   or the weighted number of records up to time \eqn{T} of the mth series,
#'   \eqn{N_{T,m}^\omega = \sum_{t=2}^{T} \omega_t I_{t,m}}; for the upper and
#'   lower records and for the forward and backward series of the
#'   \eqn{(X_{t,1}),\ldots,(X_{t,M})}, \eqn{t=1,\ldots,T}, series are used for
#'   the statistics:
#'   
#'   If \code{statistic == "d", distribution = "normal"}, 
#'   \deqn{dM^\omega = \sum_{m=1}^{M} \sum_{t=2}^{T} \omega_t\,\Big( I_{t,m} - I_{t,m}^L\Big).}
#'
#'   If \code{statistic == "D", distribution = "normal"}, 
#'   \deqn{DM^\omega = \sum_{m=1}^{M} \sum_{t=2}^{T} \Big( I_{t,m} - I_{t,m}^L - I_{t,m}^{rev} + I_{t,m}^{L,rev}\Big).}
#'
#'   If \code{statistic == "TM", distribution = "normal"}, 
#'   \deqn{TM^\omega = \sum_{m=1}^{M} \sum_{t=1}^{T} \Big( I_{t,m} + I_{t,m}^{L,rev}\Big).}
#'  
#'   While their means are very simple to calculate, their variances become 
#'   unwieldly expressions and are given by Diersen and Trenkler (2001). 
#'   
#'   The p-value is calculated with the Normal asymptotic distribution in the 
#'   usual way with a one-sided alternative based on the \code{trend} value.
#'   
#'   If \code{distribution = "t"}, for the above statistics it is calculated a
#'   new one, e.g. for \eqn{DM^\omega},
#'   \deqn{DM_S^\omega = \frac{DM^\omega - E(DM^\omega)}{\sqrt{\widehat{Var}(DM^\omega)}},}
#'   where \eqn{E(DM^\omega)} is the expectation under the null, and 
#'   \eqn{\widehat{Var}(DM^\omega)} is an estimation of the variance obtained 
#'   from the sample, which is why \eqn{M>1} is required. The statistic above
#'   is asymptotically t distributed and it has been proved that 
#'   \eqn{DM_S^\omega} is highly robust against serial correlation.
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param statistic A character string indicating the type of statistic to be 
#'   calculated, i.e., \code{"D"} (for \eqn{DM}, \eqn{DM^\omega} or 
#'   \eqn{DM^\omega_S}), \code{"d"} (for \eqn{dM}, \eqn{dM^\omega} or 
#'   \eqn{dM^\omega_S}) and \code{"TM"} (for \eqn{TM}, \eqn{TM^\omega} or 
#'   \eqn{TM^\omega_S}).
#' @param distribution A character string indicating the asymptotic 
#'   distribution of the statistic, Normal distribution "normal" or 
#'   Student's t-distribution "t".
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param trend A character string indicating the type of alternative 
#'   hypothesis, positive trend in location \code{"positive"} or
#'   negative trend in location \code{"negative"}.
#' @return A \code{"htest"} object with elements:
#'   \item{statistic}{Value of the test statistic.}
#'   \item{parameter}{Only if \code{distribution = 't'}, degrees of freedom of 
#'     \eqn{t} statistic equal to \eqn{M-1}.}
#'   \item{p.value}{P-value.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{data.name}{A character string giving the name of the data.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{foster.plot}}
#' @references
#' Diersen J, Trenkler G (1996). Records Tests for Trend in Location.
#' \emph{Statistics}, \strong{28}(1), 1-12.
#' doi: \href{https://doi.org/10.1080/02331889708802543}{10.1080/02331889708802543}
#' 
#' Diersen J, Trenkler G (2001). 
#' \emph{Weighted record tests for splitted series of observations}. 
#' En J Kunert, G Trenkler (Eds.), 
#' Mathematical Statistics with Applications in Biometry. 
#' Festschrift in Honour of Prof. Dr. Siegfried Schach (pp. 163-178). 
#' Lohmar: Josef Eul Verlang. 
#' 
#' Foster FG, Stuart A (1954). 
#' Distribution-Free Tests in Time-Series Based on the Breaking of Records.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 
#' \strong{16}(1), 1-22.
#' @examples
#' foster.test(ZaragozaSeries)
#' foster.test(ZaragozaSeries, statistic = 'd')
#' foster.test(ZaragozaSeries, weights = function(t) t-1)
#' foster.test(ZaragozaSeries, distribution = 't')
#' foster.test(ZaragozaSeries, distribution = 't', weights = function(t) t-1)
#' 
#' @export foster.test

foster.test <- function(XM_T, 
                        statistic = c('D', 'd', 'TM'), 
                        distribution = c('normal', 't'), 
                        weights = function(t) 1,
                        trend = c('positive', 'negative')) {
 
  if (class(weights) != 'function') stop("weights has to be a function")
  statistic    <- match.arg(statistic)
  distribution <- match.arg(distribution)
  trend        <- match.arg(trend)
  
  alternative <- trend == 'negative'
  
  DNAME <- deparse(substitute(XM_T))
  fun   <- deparse(weights)[2]
  
  XM_T  <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  Mcols <- ncol(XM_T) 
  w <- weights(1:Trows)
  
  if (statistic == 'D') {
    
    XM_Trev <- apply(XM_T, 2, rev)
    
    I  <- I.record.matrix(XM_T, record = 'upper') #upper records
    IL <- I.record.matrix(XM_T, record = 'lower') #lower records
    Irev  <- I.record.matrix(XM_Trev, record = 'upper') #upper records backwards
    ILrev <- I.record.matrix(XM_Trev, record = 'lower') #lower records backwards
    
    I <- I - IL - Irev + ILrev
    
    if (distribution == 'normal') {
      
      if (all(w == 1)) {
        
        DM <- sum(I)
        S  <- sqrt(4 * Mcols * sum(((3:(Trows+1) + (1:(Trows-1)) / choose(Trows, 2:Trows)) / (2:Trows)^2)))
        
        METHOD <- paste("Foster-Stuart D test for", trend, "trend in location")
        names(DM) <- "DM"
        
      } else {
        
        DM <- sum(rowSums(I) * weights(1:Trows))
        S  <- sum(weights(2:Trows)^2 / (2:Trows)) +
              sum(weights(2:(Trows - 1)) * weights((Trows - 1):2) * (exp(lfactorial(1:(Trows - 2)) + lfactorial((Trows-2):1) - lfactorial(Trows)) - 1 / Trows))
              
        for (t in 1:(Trows - 1)) {
          for (s in (t + 1):Trows) {
            
              S <- S + weights(1:Trows)[s] * weights(1:Trows)[Trows - t + 1] *
                sum(exp(lfactorial(t - 1) + lfactorial(Trows - 1:t) - log(s - 1:t) - lfactorial(t - 1:t) - lfactorial(Trows)))
          }
        }
        
        S <- sqrt(4 * Mcols * S)
        
        METHOD <- paste("Diersen-Trenkler D test for", trend, "trend in location, weights =", fun)
        names(DM) <- "DM^w"
      }
      
      pv <- stats::pnorm(q = DM, sd = S, lower.tail = alternative)
      
      structure(list(statistic = DM, p.value = pv, method = METHOD, 
                     data.name = DNAME), class = 'htest')
      
    } else { # distribution == 't'
      
      DM <- colSums(sweep(I, MARGIN = 1, weights(1:Trows), `*`))
      
      DM <- mean(DM) * sqrt(Mcols) / stats::sd(DM)
      
      Mcols <- Mcols - 1
      
      pv <- stats::pt(q = DM, df = Mcols, lower.tail = alternative)
      
      METHOD <- paste("t D test for", trend, "trend in location, weights =", fun)
      names(DM) <- "t"
      names(Mcols) <- "df"
      
      structure(list(statistic = DM, parameter = Mcols, p.value = pv, 
                     method = METHOD, data.name = DNAME), class = 'htest')
    }
  } else if (statistic == 'd') {
    
    I  <- I.record.matrix(XM_T, record = 'upper') #upper records
    IL <- I.record.matrix(XM_T, record = 'lower') #lower records
    
    I <- I - IL
    
    if (distribution == 'normal') {
      
      dM <- sum(rowSums(I) * weights(1:Trows))
      
      S <- sqrt(2 * Mcols * sum(weights(2:Trows)^2 / (2:Trows)))
      
      pv <- stats::pnorm(q = dM, sd = S, lower.tail = alternative)
      
      if (all(w == 1)) {
        
        METHOD <- paste("Foster-Stuart d test for", trend, "trend in location")
        names(dM) <- "dM"
        
      } else {
        
        METHOD <- paste("Diersen-Trenkler d test for ", trend, "trend in location, weights =", fun)
        names(dM) <- "dM^w"
      }
      
      structure(list(statistic = dM, p.value = pv, method = METHOD, 
                     data.name = DNAME), class = 'htest')
      
    } else { # distribution == 't'
      
      dM <- colSums(sweep(I, MARGIN = 1, weights(1:Trows), `*`))
      
      dM <- mean(dM)*sqrt(Mcols) / stats::sd(dM)
      
      Mcols <- Mcols - 1
      
      pv <- stats::pt(q = dM, df = Mcols, lower.tail = alternative)
      
      METHOD <- paste("t d test for", trend, "trend in location, weights =", fun)
      names(dM) <- "t"
      names(Mcols) <- "df"
      
      structure(list(statistic = dM, parameter = Mcols, p.value = pv, 
                     method = METHOD, data.name = DNAME), class = 'htest')
    }
  } else { # statistic == 'TM'
    
    I  <- I.record.matrix(XM_T, record = 'upper') #upper records
    ILrev <- I.record.matrix(apply(XM_T, 2, rev), record = 'lower') #lower records backwards
    
    I <- I + ILrev
    
    if (distribution == 'normal') {
      
      if (all(w == 1)) {
        
        TM <- sum(I)
        average <- 2 * Mcols * sum(1 / 1:Trows)
        S  <- sum(1 / 2:Trows - 1 / (2:Trows)^2) +
          sum(exp(lfactorial(1:(Trows - 2)) + lfactorial((Trows-2):1) - lfactorial(Trows)) - 1 / ((2:(Trows - 1)) * ((Trows - 1):2)))
        
        for (t in 1:(Trows - 1)) {
          for (s in (t + 1):Trows) {
            
            S <- S + sum(exp(lfactorial(t - 1) + lfactorial(Trows - 1:t) - log(s - 1:t) - lfactorial(t - 1:t) - lfactorial(Trows))) - 1 / (s * (Trows - t + 1))
          }
        }
        
        S <- sqrt(2 * Mcols * S)
        
        METHOD <- paste("Diersen-Trenkler TM test for", trend, "trend in location")
        names(TM) <- "TM"
        
      } else {
        
        TM <- sum(rowSums(I) * weights(1:Trows))
        average <- 2 * Mcols * sum(weights(1:Trows) / 1:Trows)
        S  <- sum(weights(2:Trows)^2 * (1 / 2:Trows - 1 / (2:Trows)^2)) +
          sum(weights(2:(Trows - 1)) * weights((Trows - 1):2) * (exp(lfactorial(1:(Trows - 2)) + lfactorial((Trows-2):1) - lfactorial(Trows)) - 1 / ((2:(Trows - 1)) * ((Trows - 1):2))))
        
        for (t in 1:(Trows - 1)) {
          for (s in (t + 1):Trows) {
            
            S <- S + weights(1:Trows)[s] * weights(1:Trows)[Trows - t + 1] *
              (sum(exp(lfactorial(t - 1) + lfactorial(Trows - 1:t) - log(s - 1:t) - lfactorial(t - 1:t) - lfactorial(Trows))) - 1 / (s * (Trows - t + 1)))
          }
        }
        
        S <- sqrt(2 * Mcols * S)
        
        METHOD <- paste("Diersen-Trenkler TM test for", trend, "trend in location, weights =", fun)
        names(TM) <- "TM^w"
      }
      
      pv <- stats::pnorm(q = TM, mean = average, sd = S, lower.tail = alternative)
      
      structure(list(statistic = TM, p.value = pv, method = METHOD, 
                     data.name = DNAME), class = 'htest')
      
    } else { # distribution == 't'
      
      TM <- colSums(sweep(I, MARGIN = 1, weights(1:Trows), `*`))
      average <- 2 * sum(weights(1:Trows) / 1:Trows)
      
      TM <- (mean(TM) - average) * sqrt(Mcols) / stats::sd(TM)
      
      Mcols <- Mcols - 1
      
      pv <- stats::pt(q = TM, df = Mcols, lower.tail = alternative)
      
      METHOD <- paste("t TM test for", trend, "trend in location, weights =", fun)
      names(TM) <- "TM"
      names(Mcols) <- "df"
      
      structure(list(statistic = TM, parameter = Mcols, p.value = pv, 
                     method = METHOD, data.name = DNAME), class = 'htest')
    }
  }
}