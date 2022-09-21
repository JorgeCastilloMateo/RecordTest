#' @title Number of Records Plot
#' @importFrom ggplot2 ggplot aes theme_bw labs geom_line theme
#'   scale_colour_manual scale_shape_manual scale_linetype_manual geom_ribbon 
#'   geom_errorbar geom_point annotate guides guide_legend
#' @importFrom stats qnorm
#' @description This function builds a ggplot object to compare the sample
#'   means of the (weighted) number of records in a vector up to time \eqn{t}, 
#'   \eqn{\bar N_{t}^\omega}, and the expected values 
#'   \eqn{\textrm{E}(N_t^{\omega})} 
#'   under the classical record model (i.e., of IID continuous RVs).
#' @details 
#'   This plot is associated to the test \code{\link{N.test}}.
#'   It calculates the sample means of the number of records in a set of
#'   vectors up to every time \eqn{t} (see \code{\link{Nmean.record}}). 
#'   These sample means \eqn{\bar N_{t}^\omega} are calculated from the sample of
#'   \eqn{M} values obtained from \eqn{M} vectors, the columns of matrix 
#'   \code{X}. Then, these values are plotted and compared with the expected 
#'   values \eqn{\textrm{E}(N_t^{\omega})} and their reference intervals (RIs), under
#'   the hypothesis of the classical record model. The RIs of 
#'   \eqn{\bar N_{t}^\omega} uses the fact that, under the classical record 
#'   model, the statistic is asymptotically Normal.
#'   
#'   The plot can show the four types of record at the same time (i.e., 
#'   forward upper, forward lower, backward upper and backward lower).
#'   In their interpretations one must be careful, for forward records 
#'   each time \eqn{t} corresponds to the same year of observation, but for 
#'   the backward series, time \eqn{t} corresponds to the year of observation
#'   \eqn{T-t+1} where \eqn{T} is the total number of observations in every 
#'   series. Two types of backward records can be considered (see argument
#'   \code{backward}).
#'   
#'   More details of this plot are shown in Cebrián, Castillo-Mateo, Asín 
#'   (2022).
#'   
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param record Logical vector. Vector with four elements indicating if 
#'   forward upper, forward lower, backward upper and backward lower are going
#'   to be shown, respectively. Logical values or 0,1 values are accepted.
#' @param backward A character string \code{"T"} or \code{"t"} indicating if 
#'   the backward number of records shown are calculated up to time \eqn{t} in 
#'   the backward series \eqn{\{X_T,\ldots,X_1\}} or in the series 
#'   \eqn{\{X_t,\ldots,X_1\}}. While the first option considers the evolution 
#'   of a series of records observed up to time \eqn{T}, the second considers 
#'   that until each time \eqn{t} the series has only been observed up to 
#'   \eqn{t}.
#' @param point.col,point.shape Vector with four elements indicating the colour
#'   and shape of the points. Every one of the four elements represents forward
#'   upper, forward lower, backward upper and backward lower, respectively.
#' @param conf.int Logical. Indicates if the RIs are also shown.
#' @param conf.level (If \code{conf.int == TRUE}) Confidence level of the RIs.
#' @param conf.aes (If \code{conf.int == TRUE}) A character string indicating 
#'   the aesthetic to display for the RIs, \code{"ribbon"} (grey area) or 
#'   \code{"errorbar"} (vertical lines).
#' @param conf.col Colour used to plot the expected value and (if 
#'   \code{conf.int == TRUE}) RIs.
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{N.record}}, \code{\link{N.test}}, 
#'   \code{\link{foster.test}}, \code{\link{foster.plot}}
#' @references 
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2): 313-330. 
#' \doi{10.1007/s00477-021-02122-w}
#' 
#' @examples
#' # Plot at Zaragoza, with linear weights and error bar as RIs aesthetic
#' N.plot(ZaragozaSeries, weights = function(t) t-1, conf.aes = "errorbar")
#' 
#' # Plot only upper records
#' N.plot(ZaragozaSeries, record = c(1, 0, 1, 0))
#' 
#' # Change point colour and shape
#' Zplot <- N.plot(ZaragozaSeries, 
#'   point.col = c("red", "red", "blue", "blue"), 
#'   point.shape = c(19, 4, 19, 4))
#' 
#' ## Not run: Load package ggplot2 to change the plot
#' #library("ggplot2")
#' ## Remove legend
#' #Zplot + ggplot2::theme(legend.position = "none")
#' ## Fancy axis
#' # Zplot + 
#' #   ggplot2::scale_x_continuous(name = "Year (forward)",
#' #     breaks = c(10, 30, 50, 70), 
#' #     labels=c("1960", "1980", "2000", "2020"), 
#' #     sec.axis = ggplot2::sec_axis(~ nrow(ZaragozaSeries) - . + 1951, name = "Year (backward)")) +
#' #   ggplot2::theme(axis.title.x = ggplot2::element_text(colour = "red"), 
#' #     axis.text.x = ggplot2::element_text(colour = "red"),
#' #     axis.title.x.top = ggplot2::element_text(colour = "blue"), 
#' #     axis.text.x.top = ggplot2::element_text(colour = "blue"))
#' @export N.plot
N.plot <- function(X, 
                   weights = function(t) 1, 
                   record = c("FU" = 1, "FL" = 1, "BU" = 1, "BL" = 1),
                   backward = c("T", "t"),
                   point.col = c("FU" = "red", "FL" = "blue", "BU" = "red", "BL" = "blue"),
                   point.shape = c("FU" = 19, "FL" = 19, "BU" = 4, "BL" = 4),
                   conf.int = TRUE,
                   conf.level = 0.9, 
                   conf.aes = c("ribbon", "errorbar"), 
                   conf.col = "grey69") {

  # Intro
  if (!is.function(weights)) { stop("'weights' should be a function") }
  backward <- match.arg(backward)
  conf.aes <- match.arg(conf.aes)
  record   <- as.logical(record)
  
  DNAME <- deparse(substitute(X))
  fun   <- deparse(weights)[2]
  
  X     <- as.matrix(X)
  Trows <- nrow(X)
  Mcols <- ncol(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  t     <- 1:Trows
  w     <- weights(t)
  
  if (all(w == 1)) {
    METHOD <- "Number of records"
    ylabel <- ifelse(Mcols == 1, "Number of records", "Mean number of records")
  } else {
    METHOD <- paste("Weighted number of records with weights =", fun)
    ylabel <- ifelse(Mcols == 1, "Weighted number of records", "Mean weighted number of records")
  }
  ###################################  
  # Statistic
  NmeanF.fun <- function(X, record) {
    return(cumsum(w * rowMeans(I.record(X, record = record))))
  } 
  
  NmeanB.fun <- function(X, record) {
    if (backward == "T") {
      return(cumsum(w * rowMeans(I.record(X, record = record))))
    } else {
      N <- rep(w[1], Trows)
      if (length(w) == 1) w <- rep(w, Trows)
      for (t in 2:Trows) {
        N[t] <- sum(w[1:t] * rowMeans(I.record(X[(Trows-t+1):Trows,], record = record)))
      }
      return(N)
    }  
  } 
  
  Xrev <- series_rev(X)
  df <- data.frame(t)

  if (record[1]) { N.FU <- NmeanF.fun(X, record = "upper");    df$N.FU <- N.FU }
  if (record[2]) { N.FL <- NmeanF.fun(X, record = "lower");    df$N.FL <- N.FL }
  if (record[3]) { N.BU <- NmeanB.fun(Xrev, record = "upper"); df$N.BU <- N.BU }
  if (record[4]) { N.BL <- NmeanB.fun(Xrev, record = "lower"); df$N.BL <- N.BL }
  ###################################  
  # Normal CI
  df$mu <- cumsum(w / t)
  
  if (conf.int) {
    mu <- CI1 <- CI2 <- NULL
    df$sigma <- sqrt(cumsum(w^2 * (t - 1) / t^2) / Mcols)
    q     <- stats::qnorm((1 - conf.level) / 2)
    
    df$CI1 <- df$mu + q * df$sigma   #lower bound
    df$CI2 <- df$mu - q * df$sigma   #upper bound
  }
  ###################################
  # ggplot2
  graf <- 
    ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = t)) +
    ggplot2::theme_bw() + 
    ggplot2::labs(title = METHOD, subtitle = paste("Data:", DNAME), x = "Time", y = ylabel) +
    ggplot2::geom_line(ggplot2::aes(y = mu, linetype = "CI"), colour = conf.col) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_colour_manual(name = "Records", values = point.col[record], label = c("FU" = "Forward Upper", "FL" = "Forward Lower", "BU" = "Backward Upper", "BL" = "Backward Lower")[record], breaks = c("FU", "FL", "BU", "BL")[record]) +
    ggplot2::scale_shape_manual(name = "Records", values = point.shape[record], label = c("FU" = "Forward Upper", "FL" = "Forward Lower", "BU" = "Backward Upper", "BL" = "Backward Lower")[record], breaks = c("FU", "FL", "BU", "BL")[record]) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 2))
  ###################################
  # plot CI
  if (conf.int && conf.aes == "ribbon") {
    graf <- graf +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI1, ymax = CI2, linetype = "CI"), alpha = 0.1, colour = conf.col) +
      ggplot2::scale_linetype_manual(name = "Null hyp. IID", values = c("CI" = 1), label = c("CI" = paste0("Expectation and ", 100 * conf.level, "% Normal RI")))
  } else if (conf.int && conf.aes == "errorbar") { 
    graf <- graf +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI1, ymax = CI2, linetype = 'CI'), width = 0.2, colour = conf.col) +
      ggplot2::scale_linetype_manual(name = "Null hyp. IID", values = c("CI" = 1), label = c("CI" = paste0("Expectation and ", 100 * conf.level, "% Normal RI")))
  } else {
    graf <- graf +
      ggplot2::scale_linetype_manual(name = "Null hyp. IID", values = c("CI" = 1), label = c("CI" = "Expectation"))
  }
  ###################################
  # plot points
  if (record[4]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = N.BL, colour = "BL", shape = "BL")) }
  if (record[3]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = N.BU, colour = "BU", shape = "BU")) }
  if (record[2]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = N.FL, colour = "FL", shape = "FL")) }
  if (record[1]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = N.FU, colour = "FU", shape = "FU")) }
  ###################################
  
  return(graf)
}
