#' @title Plot of the mean number of records up to time t
#' @importFrom ggplot2 ggplot aes theme_bw labs geom_line theme
#'   scale_colour_manual scale_linetype_manual geom_ribbon geom_errorbar
#'   geom_point
#' @importFrom stats qnorm
#' @description This function constructs a ggplot object to compare the sample
#'   means of the number of records in a vector up to time \eqn{t}, 
#'   \eqn{\bar N_t}, and the expected values \eqn{E(N_t)} under the classical 
#'   record model.
#' @details 
#'   First, this function calculates the sample means of the number of records
#'   in a vector up to time \eqn{t}. These sample means \eqn{\bar N_t} are 
#'   calculated from the sample of \eqn{M} values obtained from \eqn{M} 
#'   vectors, the columns of matrix \code{XM_T}. Then, these values are plotted
#'   and compared with the expected values \eqn{E(N_t)} and their confidence 
#'   intervals (CI), under the hypothesis of the classical record model.
#'
#'   The CI of \eqn{E(N_t)} uses the fact that, under the classical record 
#'   model, the statistic \eqn{\bar N_t} is asymptotically Normal.
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param record A character string indicating the type of record, 
#'   \code{"upper"}, \code{"lower"} or \code{"both"}.
#' @param interval A character string indicating the type of display of 
#'   the confidence intervals, \code{"ribbon"} (grey area) or 
#'   \code{"errorbar"} (vertical lines).
#' @param conf Confidence level of the confidence intervals.
#' @param colour Colour used to plot the expected values of \eqn{N_t}, 
#'   and the CI.
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{N_normal.test}}, \code{\link{foster.test}}
#' @examples
#' N.plot(ZaragozaSeries, weights = function(t) t-1, interval = 'errorbar')
#' 
#' Zplot <- N.plot(ZaragozaSeries, record = 'lower')
#' Zplot
#' 
#' library(ggplot2)
#' # change the colour of the points
#' colour_point <- 'black'
#' Zplot + scale_colour_manual(name = '', values = c('Lower record' = colour_point)) 
#' 
#' # remove legend
#' Zplot + theme(legend.position = 'none')
#'
#' @export N.plot

N.plot <- function(XM_T, weights = function(t) 1, 
                   record = c('both', 'upper', 'lower'), 
                   interval = c('ribbon', 'errorbar'), conf = 0.95, 
                   colour = 'salmon') {

  if (class(weights) != 'function') stop("weights has to be a function")
  record   <- match.arg(record)
  interval <- match.arg(interval)
  
  DNAME <- deparse(substitute(XM_T))
  fun   <- deparse(weights)[2]
  
  XM_T  <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  Mcols <- ncol(XM_T)
  t <- 1:Trows
  w <- weights(t)
  
  if (all(w == 1)) {
    
    METHOD <- "Number of records"
    ylabel <- "Mean number of records"
    
  } else {
    
    METHOD <- paste("Weighted number of records, weights =", fun)
    ylabel <- "Mean number of weighted records"
  }
  
  # Normal CI
  expectation <- cumsum(w / t)
  variance <- cumsum(w^2 / t * (1 - 1 / t))

  sigma <- sqrt(variance / Mcols)
    
  CI <- matrix(nrow = 2, ncol = Trows)
    
  q <- stats::qnorm((1 - conf) / 2)
    
  CI[1,] <- expectation + q * sigma   #lower bound
  CI[2,] <- expectation - q * sigma   #upper bound
  ###################################

  Nmean.fun <- function(XM_T, record) {
    I  <- I.record(XM_T, record = record)
    return(cumsum(rowMeans(I) * w))
  } 
  
  # ggplot2
  if (record=='both') {
    
    N_upper <- Nmean.fun(XM_T, record = 'upper')
    N_lower <- Nmean.fun(XM_T, record = 'lower')

    graf <- 
      ggplot2::ggplot(data = data.frame(N_upper, N_lower), mapping = ggplot2::aes(x = t)) +
      ggplot2::theme_bw() + 
      ggplot2::labs(title = METHOD, subtitle = paste('Data:', DNAME), x = "Times", y = ylabel) +
      ggplot2::geom_line(ggplot2::aes(y = expectation, linetype = 'CI'), colour = colour) +
      ggplot2::theme(legend.position = 'bottom') +
      ggplot2::scale_colour_manual(name = "", values = c('Upper record' = 'black', 'Lower record' = 'red')) +
      ggplot2::scale_linetype_manual(name = "", values = c('CI' = 1))
  
  } else {
    
    N <- Nmean.fun(XM_T, record = record)

    graf <- 
      ggplot2::ggplot(data = data.frame(N), ggplot2::aes(x = t, y = N)) + 
      ggplot2::theme_bw() + 
      ggplot2::labs(title = METHOD, subtitle = paste('Data:', DNAME), x = "Time", y = ylabel) + 
      ggplot2::geom_line(ggplot2::aes(x = t, y = expectation, linetype = 'CI'), colour = colour) +
      ggplot2::theme(legend.position = 'bottom') +
      ggplot2::scale_colour_manual(name = "", values = c('Upper record' = 'black', 'Lower record' = 'red')) +
      ggplot2::scale_linetype_manual(name = "", values = c('CI' = 1))
  }
  ###################################
  
  # plot CI
  if (interval=='ribbon') {
    
    graf <- graf +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), alpha = 0.05, colour = colour)
    
  } else { # interval=='errorbar'
    
    graf <- graf +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), width = 0.2, colour = colour)
  }
  ###################################
  
  # plot points
  if (record=='both') {
    
    graf <- graf + 
      ggplot2::geom_point(ggplot2::aes(y = N_lower, colour = 'Lower record')) +
      ggplot2::geom_point(ggplot2::aes(y = N_upper, colour = 'Upper record'))
    
  } else {
    
    colour_point <- paste(toupper(substr(record, 1, 1)), substr(record, 2, nchar(record)), ' record', sep="")
    
    graf <- graf + 
      ggplot2::geom_point(ggplot2::aes(colour = colour_point))
  } 
  ###################################
  
  return(graf)
}
