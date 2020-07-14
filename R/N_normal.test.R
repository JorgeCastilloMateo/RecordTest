#' @title Normal test on the number of records
#' @importFrom stats pnorm pt
#' @description This function performs a test based on \eqn{N_T}, the number 
#'   of records in the observation period, to study the hypothesis of the 
#'   classical record model.
#' @details 
#'   In this test the null hypothesis is that the expected value of 
#'   \eqn{N_T}, the number of records in the observation period 
#'   \eqn{t=1,\ldots,T} is \eqn{\mu_T =\sum_{i=1}^T 1/i} and the variance 
#'   \eqn{\sigma^2_T=\sum_{i=1}^T (1/i-1/i^2)}; these are the values obtained 
#'   under the classical record model (see  Arnold et al. (1998)).
#'   The test statistic is based on \eqn{\bar N_T}, the mean of the number of 
#'   records up to time t, calculated from a sample of \eqn{M} vectors 
#'   (columns in \code{XM_T}). The distribution of the statistic under the 
#'   null is asymptotically Normal.
#'
#'   If the sequences of variables in the vectors are not i.i.d., but they have
#'   a monotonous positive trend, an unilateral alternative hypothesis must 
#'   be stated, which in the case of upper records is 
#'   \eqn{\mu_T >\sum_{i=1}^T 1/i} and 
#'   \eqn{\sigma^2_T>\sum_{i=1}^T (1/i-1/i^2)}, and in the case of lower 
#'   records is \eqn{\mu_T <\sum_{i=1}^T 1/i} and 
#'   \eqn{\sigma^2_T<\sum_{i=1}^T (1/i-1/i^2)}. The opposite happens if 
#'   the trend is negative.
#'   
#'   The statistic is calculated as
#'   \deqn{Z = \frac{\bar N_T - E(\bar N_T)}{\sqrt{Var(\bar N_T)}},}
#'   which follows a standard Normal distribution.
#'   
#'   If \code{distribution = "t"}, an estimation of the variance
#'   \eqn{\widehat{Var}(DM^\omega)} is computed from the sample, which is 
#'   why \eqn{M>1} is required in this case, and the statistic is 
#'   asymptotically t distributed.
#'   
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param record A character string indicating the type of record to be 
#'   calculated, "upper" or "lower".
#' @param distribution A character string indicating the asymptotic 
#'   distribution of the statistic, Normal distribution "normal" or 
#'   Student's t-distribution "t".
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
#' @seealso \code{\link{N.plot}}
#' @references
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. New York: Wiley.
#' @examples
#' N_normal.test(ZaragozaSeries)
#' @export N_normal.test

N_normal.test <- function(XM_T, weights = function(t) 1, 
                          record = c('upper', 'lower'), 
                          distribution = c('normal', 't'),
                          trend = c('positive', 'negative')) {
  
  if(class(weights) != 'function') stop("weights has to be a function")
  record       <- match.arg(record)
  distribution <- match.arg(distribution)
  trend        <- match.arg(trend)
      
  alternative <- trend == 'negative'
  
  DNAME <- deparse(substitute(XM_T))

  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  t <- 1:Trows
  w <- weights(t)

  NT <- colSums(sweep(I.record(XM_T, record = record), MARGIN = 1, STATS = w, FUN = `*`))
  mu <- sum(w / t)
  sigma <- ifelse(distribution == 'normal', sqrt(sum(w^2 / t * (1 - 1 / t))), sd(NT))
  NT <- (mean(NT) - mu) / (sigma / sqrt(Mcols))
  
  if (distribution == 'normal') {
    if (record=='lower') {
      METHOD <- paste("Number of records normality test. Lower records and", trend, "trend")
      pvalue <- stats::pnorm(NT, mean = 0, sd = 1, lower.tail = !alternative)
    } else { # record=='upper'
      METHOD <- paste("Number of records normality test. Upper records and", trend, "trend")
      pvalue <- stats::pnorm(NT, mean = 0, sd = 1, lower.tail = alternative)
    }
    
    names(NT) <- "Z"
    
    structure(list(statistic = NT, p.value = pvalue, 
                   method = METHOD, data.name = DNAME), class='htest')
    
  } else { # distribution == 't'
    
    Mcols <- Mcols - 1
    
    if (record=='lower') {
      METHOD <- paste("Number of records t test. Lower records and", trend, "trend")
      pvalue <- stats::pt(NT, df = Mcols, lower.tail = !alternative)
    } else { # record=='upper'
      METHOD <- paste("Number of records t test. Upper records and", trend, "trend")
      pvalue <- stats::pt(NT, df = Mcols, lower.tail = alternative)
    }
    
    names(NT) <- "t"
    names(Mcols) <- "df"
    
    structure(list(statistic = NT, parameter = Mcols, p.value = pvalue, 
                   method = METHOD, data.name = DNAME), class='htest')
  }
}
