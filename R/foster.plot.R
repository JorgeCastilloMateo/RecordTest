#' @title Plot based on Foster-Stuart and Diersen-Trenkler statistics
#' @importFrom ggplot2 ggplot aes geom_point theme_bw geom_line theme
#'   geom_ribbon geom_errorbar labs
#' @importFrom stats qnorm
#' @description This function constructs a ggplot object to display two-sided
#'   confidence intervals based on Foster-Stuart and Diersen-Trenkler statistics
#'   for randomness.
#'
#' @details See \code{\link{foster.test}}.
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param statistic A character string indicating the statistic to be 
#'   calculated, i.e., \code{"D"} (for \eqn{\overline{DM}} or 
#'   \eqn{\overline{DM}^\omega}), \code{"d"} (for \eqn{\overline{dM}} or 
#'   \eqn{\overline{dM}^\omega}) and \code{"TM"} (for \eqn{\overline{TM}} or 
#'   \eqn{\overline{TM}^\omega}).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param interval A character string indicating the type of display of the 
#'   confidence intervals, \code{"ribbon"} (grey area) or \code{"errorbar"} 
#'   (vertical lines).
#' @param conf Numeric value in \eqn{(0,1)}. Confidence level of the two-sided
#'   confidence intervals.
#' @param colour Colour used to plot the expected values and the CI.
#' @return A ggplot graph object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{foster.test}}, \code{\link{N.plot}}, 
#'   \code{\link{N_normal.test}}
#' @references
#' Diersen J, Trenkler G (1996). Records Tests for Trend in Location.
#' \emph{Statistics}, \strong{28}(1), 1-12.
#' doi:\href{https://doi.org/10.1080/02331889708802543}{10.1080/02331889708802543}
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
#' foster.plot(ZaragozaSeries)
#' foster.plot(ZaragozaSeries, interval = 'error', conf = 0.9, colour = 1)
#' foster.plot(ZaragozaSeries, statistic = 'd', weights = function(t) t-1)
#' 
#' @export foster.plot

foster.plot <- function(XM_T, statistic = c('D', 'd', 'TM'), 
                        weights = function(t) 1,
                        interval = c('ribbon', 'errorbar'), 
                        conf = 0.95, colour = 'salmon') {
  
  if (class(weights) != 'function') stop("weights has to be a function")
  statistic <- match.arg(statistic)
  interval  <- match.arg(interval)
  
  DNAME <- deparse(substitute(XM_T))
  fun   <- deparse(weights)[2]
  
  XM_T  <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  Mcols <- ncol(XM_T) 
  w <- weights(1:Trows)
  
  if (statistic == 'D') {
    
    XM_Trev <- apply(XM_T, 2, rev)
    
    I  <- I.record(XM_T, record = 'upper') #upper records
    IL <- I.record(XM_T, record = 'lower') #lower records
    Irev  <- I.record(XM_Trev, record = 'upper')
    ILrev <- I.record(XM_Trev, record = 'lower')
    
    I <- I - IL - Irev + ILrev
    
    DM <- cumsum(rowMeans(I) * weights(1:Trows))
    
    if (all(w == 1)) {
      
      METHOD <- "Plot based on Foster-Stuart D statistic"
      ylabel <- expression(bar(N)[t] - bar(N)[t]^L - bar(N)[t]^rev + bar(N)[t]^list(L,rev))
      
      sigma <- c(0, sqrt(4 / Mcols * cumsum(((3:(Trows+1) + (1:(Trows-1)) / choose(Trows, 2:Trows)) / (2:Trows)^2))))
      
    } else {
      
      METHOD <- paste("Plot based on Diersen-Trenkler D statistic, weights =", fun)
      ylabel <- expression(bar(N)[t]^omega - bar(N)[t]^list(L,omega) - bar(N)[t]^list(rev,omega) + bar(N)[t]^list(L,rev,omega))
      
      sigma <- sqrt(c(0, sapply(X = 2:Trows, FUN = variance.fun, weights = weights, simplify = TRUE)) / Mcols)
    }
    
    average <- 0
         
  } else if (statistic == 'd') {
    
    I  <- I.record(XM_T, record = 'upper') #upper records
    IL <- I.record(XM_T, record = 'lower') #lower records
    
    I <- I - IL
    
    DM <- cumsum(rowMeans(I) * weights(1:Trows))
    
    if(all(w == 1)){
      
      METHOD <- "Plot based on Foster-Stuart d statistic"
      ylabel <- expression(bar(N)[t] - bar(N)[t]^L)
      
    } else {
      
      METHOD <- paste("Plot based on Diersen-Trenkler d statistic, weights =", fun)
      ylabel <- expression(bar(N)[t]^omega - bar(N)[t]^list(L,omega))
    }
    
    average <- 0
    sigma <- c(0, sqrt(2 / Mcols * cumsum(weights(2:Trows)^2 / (2:Trows))))
    
  } else { # statistic == 'TM'
    
    XM_Trev <- apply(XM_T, 2, rev)
    
    I  <- I.record(XM_T, record = 'upper') #upper records
    ILrev <- I.record(XM_Trev, record = 'lower')
    
    I <- I + ILrev
    
    DM <- cumsum(rowMeans(I) * weights(1:Trows))
    
    if (all(w == 1)) {
      
      METHOD <- "Plot based on Diersen-Trenkler TM statistic"
      ylabel <- expression(bar(N)[t] + bar(N)[t]^list(L,rev))
      
    } else {
      
      METHOD <- paste("Plot based on Diersen-Trenkler TM statistic, weights =", fun)
      ylabel <- expression(bar(N)[t]^omega + bar(N)[t]^list(L,rev,omega))
    }
    
    average <- 2 * cumsum(weights(1:Trows) / 1:Trows)
    sigma <- sqrt(c(0, sapply(X = 2:Trows, FUN = variance.fun, weights = weights, simplify = TRUE, statistic = statistic)) / Mcols)
  }
  
  CI <- matrix(nrow = 2, ncol = Trows)
  
  q <- stats::qnorm((1 - conf) / 2)
  
  CI[1,] <- average + q * sigma   #lower bound
  CI[2,] <- average - q * sigma   #upper bound
  
  # plot
  graf <- 
    ggplot2::ggplot(data = data.frame(DM), ggplot2::aes(x = 1:Trows, y = DM)) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::labs(title = METHOD, 
                  subtitle = paste('Data:', DNAME), 
                  x = "t", y = ylabel)
  
  if (statistic == 'TM') graf <- graf + ggplot2::geom_line(ggplot2::aes(x = 1:Trows, y = average), colour = colour)
  else graf <- graf + ggplot2::geom_line(ggplot2::aes(y = 0), colour = colour)
  ###################################
  
  # plot CI
  if (interval == 'ribbon') {
    
    graf <- graf +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), alpha = 0.05, colour = colour)
    
  } else {
    
    graf <- graf +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), width = 0.2, colour = colour)
  }
  ###################################
  
  graf <- graf + ggplot2::geom_point()
  
  return(graf)
}

variance.fun <- function(iter, weights, statistic = 'D') { # necessary sqrt(variance.fun(...) */ Mcols)
  
  Trows <- iter
  
  if (statistic == 'D') {
    
    if (Trows == 2) {
    
      S <- weights(2)^2
    
    } else { # Trows > 2
    
      S <- sum(weights(2:Trows)^2 / (2:Trows)) +
        sum(weights(2:(Trows - 1)) * weights((Trows - 1):2) * (exp(lfactorial(1:(Trows - 2)) + lfactorial((Trows-2):1) - lfactorial(Trows)) - 1 / Trows))
  
      for (t in 1:(Trows - 1)) {
        for (s in (t + 1):Trows) {
      
          S <- S + weights(1:Trows)[s] * weights(1:Trows)[Trows - t + 1] *
            sum(exp(lfactorial(t - 1) + lfactorial(Trows - 1:t) - log(s - 1:t) - lfactorial(t - 1:t) - lfactorial(Trows)))
        }
      }
    }
    S <- 4 * S
    
  } else { # statistic == 'TM'
    
    if (Trows == 2) {
      
      S <- weights(2)^2 / 2
      
    } else { # Trows > 2
      if (length(weights(1:Trows)) == 1) weights <- function(t) rep(1, length(t))
      
      S  <- sum(weights(2:Trows)^2 * (1 / 2:Trows - 1 / (2:Trows)^2)) +
        sum(weights(2:(Trows - 1)) * weights((Trows - 1):2) * (exp(lfactorial(1:(Trows - 2)) + lfactorial((Trows-2):1) - lfactorial(Trows)) - 1 / ((2:(Trows - 1)) * ((Trows - 1):2))))
    
      for (t in 1:(Trows - 1)) {
        for (s in (t + 1):Trows) {
        
          S <- S + weights(1:Trows)[s] * weights(1:Trows)[Trows - t + 1] *
            (sum(exp(lfactorial(t - 1) + lfactorial(Trows - 1:t) - log(s - 1:t) - lfactorial(t - 1:t) - lfactorial(Trows))) - 1 / (s * (Trows - t + 1)))
        }
      }
    
      S <- 2 * S
    }
  }

  return(S)
}
