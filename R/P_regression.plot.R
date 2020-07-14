#' @title Plot of record probabilities
#' @importFrom ggplot2 ggplot aes theme_bw geom_line geom_smooth labs 
#'   geom_ribbon geom_errorbar geom_point
#' @importFrom stats qbinom
#' @description This function constructs a ggplot object to display different 
#'   functions of the record probabilities at time \eqn{t}, \eqn{p_t}.
#' @details Three different types of plots which aim to analyse the hypothesis
#'   of the classical record model using the record probabilities are 
#'   implemented. Estimations of the record probabilities \eqn{\hat p_t} used
#'   in the plots are obtained as the proportion of records at time \eqn{t}
#'   in \eqn{M} vectors (columns of matrix \code{XM_T}).
#'
#'   Type 2 is the plot of the estimated record probabilities \eqn{p_t} versus
#'   time. The expected probabilities under the classical record model, 
#'   \eqn{p_t=1/t}, are also plotted, together with binomial confidence 
#'   intervals. Type 3 is the same plot but on a logarithmic scale, so that the
#'   expected value is \eqn{-log(t)}.
#'
#'   Type 1 is the plot of the observed values \eqn{t \hat p_t} versus time.
#'   The expected values under the classical record model are \eqn{1} for any
#'   value \eqn{t}, so that a cloud of points around \eqn{1} and with no trend
#'   should be expected. The estimated values are plotted, together with 
#'   binomial confidence intervals. In addition, a regression line is fitted
#'   to the cloud of points and plotted together with confidence intervals of
#'   the response. If the classical record model is true, the confidence band 
#'   (in grey) should contain the horizontal line equal to \eqn{1}. Plots of
#'   type 1 are easier to interpret than types 2 and 3.
#'
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be 
#'   calculated, \code{"upper"} or \code{"lower"}.
#' @param plot One of the values 1, 2 or 3. It determines the type of plot to
#'   be displayed. See Details.
#' @param interval A character string indicating the type of display of the
#'   confidence intervals, \code{"ribbon"} (grey area) or \code{"errorbar"}
#'   (vertical lines).
#' @param conf Numeric value in (0,1). Confidence level of the confidence
#'   intervals.
#' @param weight logical. If \code{TRUE} (the default) the regression line is
#'   estimated with weights. Only used if \code{plot = 1}.
#' @param colour_point Colour used to plot the points.
#' @param colour_CI Colour used to plot the expected values and the CI.
#' @param colour_lm Colour used to plot the regression line. Only used if
#'   \code{plot = 1}.
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{P_regression.test}}
#' @examples
#' P_regression.plot(ZaragozaSeries, plot = 1)
#' P_regression.plot(ZaragozaSeries, plot = 2)
#' P_regression.plot(ZaragozaSeries, plot = 3)
#' P_regression.plot(series_rev(ZaragozaSeries), record = 'lower')
#'
#' @export

P_regression.plot <- function(XM_T, 
                              record = c('upper', 'lower'), 
                              plot = 1,
                              interval = c('ribbon', 'errorbar'),
                              conf = 0.95, 
                              weight = TRUE, 
                              colour_point = 'black', 
                              colour_CI = 'salmon', 
                              colour_lm = 'royalblue4'){
  
  record   <- match.arg(record)
  interval <- match.arg(interval)
  
  DNAME <- deparse(substitute(XM_T))
  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  t <- 1:Trows
  
  p <- P.record(XM_T, record = record)
  
  # CI
  CI <- matrix(nrow=2, ncol=Trows)
  
  CI[1,] <- stats::qbinom((1 - conf) / 2, size = Mcols, prob = 1 / (1:Trows))
  CI[2,] <- stats::qbinom((1 + conf) / 2, size = Mcols, prob = 1 / (1:Trows))
  ###################################
  ### type 1
  if (plot == 1) { 
  CI[1,] <- t / Mcols * CI[1,]
  CI[2,] <- t / Mcols * CI[2,]
  tp <- t * p
  weights <- ifelse(rep(weight, Trows), c(0, Mcols / t[-Trows]), c(0, rep(1, Trows - 1)))

  graf <- ggplot2::ggplot(data = data.frame(t, tp), ggplot2::aes(x = t, y = tp)) +
    ggplot2::theme_bw() + 
    ggplot2::geom_line(ggplot2::aes(x = t, y = rep(1, Trows)), colour = colour_CI) +
    ggplot2::geom_smooth(formula = y ~ x, method = lm, mapping = ggplot2::aes(weight = weights), colour = colour_lm) +
    ggplot2::labs(title = paste("Regression plot on", record, "record probabilities"), subtitle = paste("Data:", DNAME), x = "t", y = expression(t %*% hat(p)[t]))
  ### type 2
  } else if (plot == 2) {
  CI[1,] <- CI[1,] / Mcols
  CI[2,] <- CI[2,] / Mcols

  graf <- ggplot2::ggplot(data = data.frame(t, p), ggplot2::aes(x = t, y = p)) +
    ggplot2::theme_bw() + 
    ggplot2::geom_line(ggplot2::aes(x = t, y = 1 / t), colour = colour_CI) +
    ggplot2::labs(title = paste("Plot based on", record, "record probabilities"), subtitle = paste("Data:", DNAME), x = "t", y = expression(hat(p)[t]))
  ### type 3
  } else if (plot == 3) {
  CI[1,] <- log(CI[1,] / Mcols)
  CI[2,] <- log(CI[2,] / Mcols)
  logp <- log(p)

  graf <- ggplot2::ggplot(data = data.frame(t, logp), ggplot2::aes(x = log(t), y = logp)) +
    ggplot2::theme_bw() + 
    ggplot2::geom_line(ggplot2::aes(x = log(t), y = -log(t)), colour = colour_CI) +
    ggplot2::labs(title = paste("Plot based on", record, "record probabilities (logarithmic scale)"), subtitle = paste("Data:", DNAME), x = "log(t)", y = expression(log(hat(p)[t])))
  } else stop("'plot' should be one of 1, 2, or 3")
  ###################################
  # adding CI
  if (interval == 'ribbon') {
    
    graf <- graf + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), alpha = 0.05, colour = colour_CI)
  
  } else {
    
    graf <- graf + 
      ggplot2::geom_errorbar(ggplot2::aes(ymin = CI[1,], ymax = CI[2,]), width = 0.2, colour = colour_CI)
  } 
  ###################################
  
  graf <- graf + ggplot2::geom_point(colour = colour_point)
  
  return(graf)
}
