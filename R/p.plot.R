#' @title Probabilities of Record Plots
#' @importFrom ggplot2 ggplot aes theme_bw geom_line geom_smooth labs 
#'   geom_ribbon geom_errorbar geom_point geom_step
#' @importFrom stats qbinom
#' @description This function builds a ggplot object to display different 
#'   functions of the record probabilities at time \eqn{t}, \eqn{p_t}.
#'   A graphical tool to study the hypothesis of the classical record model 
#'   (i.e., of IID continuous RVs).
#' @details 
#'   Three different types of plots which aim to analyse the hypothesis
#'   of the classical record model using the record probabilities are 
#'   implemented. Estimations of the record probabilities \eqn{\hat p_t} used
#'   in the plots are obtained as the proportion of records at time \eqn{t}
#'   in \eqn{M} vectors (columns of matrix \code{X}) (see 
#'   \code{\link{p.record}}).
#'
#'   Type 1 is the plot of the observed values \eqn{t \hat p_t} versus time 
#'   \eqn{t} (see \code{\link{p.regression.test}} for its associated test and
#'   details). 
#'   The expected values under the classical record model are \eqn{1} for any
#'   value \eqn{t}, so that a cloud of points around \eqn{1} and with no trend
#'   should be expected. The estimated values are plotted, together with 
#'   binomial confidence intervals (CIs). In addition, a smoothing function
#'   can be fitted to the cloud of points.
#'   
#'   Type 2 is the plot of the estimated record probabilities \eqn{p_t} versus
#'   time \eqn{t}. The expected probabilities under the classical record model, 
#'   \eqn{p_t=1/t}, are also plotted, together with binomial CIs. 
#'   
#'   Type 3 is the same plot but on a logarithmic scale, so that the
#'   expected value is \eqn{-\log(t)}. In this case, another smoothing 
#'   function can be fitted to the cloud of points.
#'   
#'   Type 1 plot was proposed by Cebrián, Castillo-Mateo, Asín (2021), while 
#'   type 2 and 3 appear in Benestad (2003, Figures 8 and 9, 2004, Figure 4).
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param plot One of the values "1", "2" or "3" (character or numeric class 
#'   are both allowed). It determines the type of plot to be displayed (see 
#'   Details).
#' @param record Logical vector. Vector with four elements indicating if 
#'   forward upper, forward lower, backward upper and backward lower are going
#'   to be shown, respectively. Logical values or 0,1 values are accepted.
#' @param point.col,point.shape Vector with four elements indicating the colour
#'   and shape of the points. Every one of the four elements represents forward
#'   upper, forward lower, backward upper and backward lower, respectively.
#' @param conf.int Logical. Indicates if the CIs are also shown.
#' @param conf.level (If \code{conf.int == TRUE}) Confidence level of the CIs.
#' @param conf.aes (If \code{conf.int == TRUE}) A character string indicating 
#'   the aesthetic to display for the CIs, \code{"ribbon"} (grey area) or 
#'   \code{"errorbar"} (vertical lines).
#' @param conf.col Colour used to plot the expected value and (if 
#'   \code{conf.int == TRUE}) CIs.
#' @param smooth (If \code{plot = 1} or \code{3}) Logical. If \code{TRUE}, a
#'   smoothing in the probabilities is also plotted.
#' @param smooth.formula (\code{smooth = TRUE}) \code{\link{formula}} to use 
#'   in the smooth function, e.g., \code{y ~ x}, 
#'   \code{y ~ poly(x, 2, raw = TRUE)}, \code{y ~ log(x)}. 
#' @param smooth.method (If \code{smooth = TRUE}) Smoothing method (function) 
#'   to use, e.g., \code{\link{lm}} or \code{\link{loess}}.
#' @param smooth.weight (If \code{smooth = TRUE}) Logical. If \code{TRUE} 
#'   (the default) the smoothing is estimated with weights.
#' @param smooth.linetype (If \code{smooth = TRUE}) Vector with four elements 
#'   indicating the line type of the smoothing. Every one of the four elements
#'   represents forward upper, forward lower, backward upper and backward
#'   lower, respectively.
#' @param ... Further arguments to pass through the smooth 
#'   (see \code{ggplot2::geom_smooth}).
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{p.regression.test}}
#' @references 
#' Benestad RE (2003). 
#' “How Often Can We Expect a Record Event?” 
#' \emph{Climate Research}, \strong{25}(1), 3–13.
#' 
#' Benestad RE (2004). 
#' “Record-Values, Nonstationarity Tests and Extreme Value Distributions.” 
#' \emph{Global and Planetary Change}, \strong{44}(1–4), 11–26. 
#' 
#' Cebrián A, Castillo-Mateo J and Asín J (2021).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' Available at Research Square \doi{10.21203/rs.3.rs-214787/v1}
#' @examples
#' # three plots available
#' p.plot(ZaragozaSeries, plot = 1)
#' p.plot(ZaragozaSeries, plot = 2)
#' p.plot(ZaragozaSeries, plot = 3)
#'
#' # Posible fits (plot 1):
#' #fit a line
#' p.plot(ZaragozaSeries, record = c(1,0,0,0))
#' # fit a second order polynomial
#' p.plot(ZaragozaSeries, record = c(1,0,0,0), 
#'   smooth.formula = y ~ poly(x, degree = 2))
#' # force the line to pass by E(t*p_t) = 1 when t = 1, i.e., E(t*p_t) = 1 + beta_1 * (t-1)
#' p.plot(ZaragozaSeries, record = c(1,0,0,0), 
#'   smooth.formula = y ~ I(x-1) - 1 + offset(rep(1, length(x))))
#' # force the second order polynomial pass by E(t*p_t) = 1 when t = 1
#' p.plot(ZaragozaSeries, record = c(1,0,0,0), 
#'   smooth.formula = y ~ I(x-1) + I(x^2-1) - 1 + offset(rep(1, length(x))))
#' # fit a loess
#' p.plot(ZaragozaSeries, record = c(1,0,0,0), 
#'   smooth.method = stats::loess, span = 0.25)
#' @export

p.plot <- function(X, 
                   plot = c("1", "2", "3"),
                   record = c("FU" = 1, "FL" = 1, "BU" = 1, "BL" = 1),
                   point.col = c("FU" = "red", "FL" = "blue", "BU" = "red", "BL" = "blue"),
                   point.shape = c("FU" = 19, "FL" = 19, "BU" = 4, "BL" = 4),
                   conf.int = TRUE,
                   conf.level = 0.9, 
                   conf.aes = c("ribbon", "errorbar"),
                   conf.col = "grey69",
                   smooth = TRUE,
                   smooth.formula = y ~ x,
                   smooth.method = stats::lm,
                   smooth.weight = TRUE, 
                   smooth.linetype = c("FU" = 1, "FL" = 1, "BU" = 2, "BL" = 2),
                   ...){
  
  # Intro
  plot <- as.character(plot)
  plot <- match.arg(plot)
  conf.aes <- match.arg(conf.aes)
  
  DNAME <- deparse(substitute(X))
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  t <- 1:Trows
  ###################################
  # Statistic
  Xrev <- series_rev(X)
  df <- data.frame(t)
  
  if (record[1]) { p.FU <- p.record(X, record = "upper");    df$p.FU <- p.FU }
  if (record[2]) { p.FL <- p.record(X, record = "lower");    df$p.FL <- p.FL }
  if (record[3]) { p.BU <- p.record(Xrev, record = "upper"); df$p.BU <- p.BU }
  if (record[4]) { p.BL <- p.record(Xrev, record = "lower"); df$p.BL <- p.BL }
  ###################################
  # Binomial CI
  if (conf.int) {
    CI1 <- CI2 <- NULL
    df$CI1 <- stats::qbinom((1 - conf.level) / 2, size = Mcols, prob = 1 / t) / Mcols
    df$CI2 <- stats::qbinom((1 + conf.level) / 2, size = Mcols, prob = 1 / t) / Mcols
  }
  ###################################
  switch(plot,
    "1" = {
        ### type 1
        if (smooth.weight) { weights <- c(0, Mcols / t[-Trows]) } 
        else               { weights <- c(0, rep(1, Trows - 1)) }
        
        df[-1] <- t * df[-1]
        df$weights <- weights
        
        graf <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = t)) +
          ggplot2::theme_bw() + 
          ggplot2::geom_line(ggplot2::aes(y = rep(1, Trows), size = "CI"), colour = conf.col) +
          ggplot2::labs(title = "Normalised probabilities of record", subtitle = paste("Data:", DNAME), x = "t", y = expression(t %*% hat(p)[t]))
    },
    "2" = {
        ### type 2
        graf <- ggplot2::ggplot(data = df, ggplot2::aes(x = t)) +
          ggplot2::theme_bw() + 
          ggplot2::geom_line(ggplot2::aes(x = t, y = 1 / t, size = "CI"), colour = conf.col) +
          ggplot2::labs(title = "Probabilities of record", subtitle = paste("Data:", DNAME), x = "t", y = expression(hat(p)[t]))
    },
    "3" = {
        ### type 3
        if (smooth.weight) {
          x <- 1:Mcols
          E <- function(t, x) { 
            return(sum(log(x / Mcols) * choose(Mcols, x) * (t - 1)^(Mcols - x) / t^Mcols) / (1 - (1 - 1/t)^Mcols))
          }
          VAR <- function(t) {
            return(sum((log(x / Mcols) - E(t, x))^2 * choose(Mcols, x) * (t - 1)^(Mcols - x) / t^Mcols) / (1 - (1 - 1/t)^Mcols))
          }
          weights <- c(0, rep(NA, Trows - 1))
          for (i in 2:Trows) { weights[i] <-  1 / VAR(i) }
        } else { weights <- c(0, rep(1, Trows - 1)) }
      
        df <- log(df)
        df$weights <- weights
        
        graf <- ggplot2::ggplot(data = df, ggplot2::aes(x = t)) +
          ggplot2::theme_bw() + 
          ggplot2::geom_line(ggplot2::aes(x = t, y = -t, size = "CI"), colour = conf.col) +
          ggplot2::labs(title = paste("Log-probabilities of record"), subtitle = paste("Data:", DNAME), x = "log(t)", y = expression(log(hat(p)[t]))) 
    }
  )
  ###################################
  # adding CI
  if (conf.int && conf.aes == "ribbon") {
    graf <- graf + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI1, ymax = CI2, size = "CI"), alpha = 0.1, colour = conf.col) +
      ggplot2::scale_size_manual(name = "Null hyp. IID", values = c("CI" = 0.5), label = c("CI" = paste0("Expectation and ", 100 * conf.level, "% Binomial CI")))
  } else if (conf.int && conf.aes == "errorbar") {
    graf <- graf + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI1, ymax = CI2, size = "CI"), width = 0.2, colour = conf.col) +
      ggplot2::scale_size_manual(name = "Null hyp. IID", values = c("CI" = 0.5), label = c("CI" = paste0("Expectation and ", 100 * conf.level, "% Binomial CI")))
  } else {
    graf <- graf +
      ggplot2::scale_size_manual(name = "Null hyp. IID", values = c("CI" = 0.5), label = c("CI" = "Expectation"))
  }
  ###################################
  # plot points
  if (record[4]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = p.BL, colour = "BL", shape = "BL")) }
  if (record[3]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = p.BU, colour = "BU", shape = "BU")) }
  if (record[2]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = p.FL, colour = "FL", shape = "FL")) }
  if (record[1]) { graf <- graf + ggplot2::geom_point(ggplot2::aes(y = p.FU, colour = "FU", shape = "FU")) }
  
  if (plot != 2 && smooth) {
    if (record[4]) { graf <- graf + ggplot2::geom_smooth(formula = smooth.formula, method = smooth.method, ..., mapping = ggplot2::aes(y = p.BL, weight = weights, colour = "BL", linetype = "BL"), se = FALSE, alpha = 0.1) }
    if (record[3]) { graf <- graf + ggplot2::geom_smooth(formula = smooth.formula, method = smooth.method, ..., mapping = ggplot2::aes(y = p.BU, weight = weights, colour = "BU", linetype = "BU"), se = FALSE, alpha = 0.1) }
    if (record[2]) { graf <- graf + ggplot2::geom_smooth(formula = smooth.formula, method = smooth.method, ..., mapping = ggplot2::aes(y = p.FL, weight = weights, colour = "FL", linetype = "FL"), se = FALSE, alpha = 0.1) }
    if (record[1]) { graf <- graf + ggplot2::geom_smooth(formula = smooth.formula, method = smooth.method, ..., mapping = ggplot2::aes(y = p.FU, weight = weights, colour = "FU", linetype = "FU"), se = FALSE, alpha = 0.1) }
  
    graf <- graf + 
      ggplot2::scale_linetype_manual(name = "Records", values = smooth.linetype, label = c("FU" = "Forward Upper", "FL" = "Forward Lower", "BU" = "Backward Upper", "BL" = "Backward Lower"), breaks = c("FU", "FL", "BU", "BL"))
  }
  ###################################
     
  graf <- graf + 
    ggplot2::scale_colour_manual(name = "Records", values = point.col, label = c("FU" = "Forward Upper", "FL" = "Forward Lower", "BU" = "Backward Upper", "BL" = "Backward Lower"), breaks = c("FU", "FL", "BU", "BL")) +
    ggplot2::scale_shape_manual(name = "Records", values = point.shape, label = c("FU" = "Forward Upper", "FL" = "Forward Lower", "BU" = "Backward Upper", "BL" = "Backward Lower"), breaks = c("FU", "FL", "BU", "BL")) +
    ggplot2::theme(legend.position = "bottom")
  
  return(graf)
}