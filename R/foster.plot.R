#' @title Plots Based on Foster-Stuart and Diersen-Trenkler Statistics
#' 
#' @importFrom ggplot2 ggplot aes geom_point theme_bw geom_line theme
#'   geom_ribbon geom_errorbar labs
#' @importFrom stats qnorm
#' 
#' @description This function builds a ggplot object to display two-sided
#'   reference intervals based on Foster-Stuart and Diersen-Trenkler 
#'   statistics to study the hypothesis of the classical record model (i.e., 
#'   of IID continuous RVs).
#'   
#' @details 
#'   The function plots the observed values of the statistic selected with 
#'   \code{statistic}, obtained with the series up to time \eqn{t} for 
#'   \eqn{t=1,\ldots,T}. The plot also includes the expected values and 
#'   reference intervals (RIs) based on the asymptotic normal distribution 
#'   of the statistics under the null hypothesis.
#' 
#'   This function implements the same ideas that \code{\link{N.plot}}, but with
#'   the statistics computed in \code{\link{foster.test}}. 
#'   
#'   These plots are useful to see the evolution in the record occurrence 
#'   and to follow the evolution of the trend. The plot was proposed by 
#'   Cebrián, Castillo-Mateo, Asín (2022) where its application is shown. 
#'   
#' @param X A numeric vector, matrix (or data frame).
#' @param weights A function indicating the weight given to the different 
#'   records according to their position in the series,
#'   e.g., if \code{function(t) t-1} then \eqn{\omega_t = t-1}.
#' @param statistic A character string indicating the type of statistic to be 
#'   calculated, i.e., one of \code{"D"}, \code{"d"}, \code{"S"}, \code{"s"},
#'   \code{"U"}, \code{"L"} or \code{"W"} (see \code{\link{foster.test}}).
#' @param point.col,point.shape Value with the colour and shape of the points. 
#' @param conf.int Logical. Indicates if the RIs are also shown.
#' @param conf.level (If \code{conf.int == TRUE}) Confidence level of the RIs.
#' @param conf.aes (If \code{conf.int == TRUE}) A character string indicating 
#'   the aesthetic to display for the RIs, \code{"ribbon"} (grey area) or 
#'   \code{"errorbar"} (vertical lines).
#' @param conf.col Colour used to plot the expected value and (if 
#'   \code{conf.int == TRUE}) RIs.
#' @return A ggplot graph object.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{foster.test}}, \code{\link{N.plot}}, 
#'   \code{\link{N.test}}
#' @references
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
#' foster.plot(ZaragozaSeries)
#' # D-statistic with linear weights
#' foster.plot(ZaragozaSeries, weights = function(t) t-1)
#' # S-statistic with linear weights
#' foster.plot(ZaragozaSeries, statistic = "S", weights = function(t) t-1)
#' # U-statistic with weights (upper tail)
#' foster.plot(ZaragozaSeries, statistic = "U", weights = function(t) t-1)
#' # L-statistic with weights (lower tail)
#' foster.plot(ZaragozaSeries, statistic = "L", weights = function(t) t-1)
#' 
#' @export foster.plot
foster.plot <- function(X, 
                        weights = function(t) 1,
                        statistic = c("D", "d", "S", "s", "U", "L", "W"),
                        point.col = "black",
                        point.shape = 19,
                        conf.int = TRUE,
                        conf.level = 0.9,
                        conf.aes = c("ribbon", "errorbar"), 
                        conf.col = "grey69") {
  
  # Intro
  if (!is.function(weights)) { stop("'weights' should be a function") }
  statistic <- match.arg(statistic)
  conf.aes  <- match.arg(conf.aes)
  
  DNAME <- deparse(substitute(X))
  X <- as.matrix(X)
  Trows <- nrow(X)
  Mcols <- ncol(X) 
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }
  t <- 1:Trows
  w <- weights(t)
  
  switch(statistic,
    "D" = {
      METHOD <- "Foster-Stuart D-statistic"
      ylabel <- ifelse(all(w == 1), 
                       "Mean D-statistic",
                       "Mean weighted D-statistic"
                       #expression(bar(N)[t]^list((FU)) - bar(N)[t]^list((FL)) - bar(N)[t]^list((BU)) + bar(N)[t]^list((BL))),
                       #expression(bar(N)[t]^list(omega,(FU)) - bar(N)[t]^list(omega,(FL)) - bar(N)[t]^list(omega,(BU)) + bar(N)[t]^list(omega,(BL)))
      )
      },
    "d" = {
      METHOD <- "Foster-Stuart d-statistic"
      ylabel <- ifelse(all(w == 1), 
                       "Mean d-statistic",
                       "Mean weighted d-statistic"
                       #expression(bar(N)[t]^list((FU)) - bar(N)[t]^list((FL))),
                       #expression(bar(N)[t]^list(omega,(FU)) - bar(N)[t]^list(omega,(FL)))
      )
      },
    "S" = {
      METHOD <- "Foster-Stuart S-statistic"
      ylabel <- ifelse(all(w == 1), 
                       "Mean S-statistic",
                       "Mean weighted S-statistic"
                       #expression(bar(N)[t]^list((FU)) + bar(N)[t]^list((FL)) + bar(N)[t]^list((BU)) + bar(N)[t]^list((BL))),
                       #expression(bar(N)[t]^list(omega,(FU)) + bar(N)[t]^list(omega,(FL)) + bar(N)[t]^list(omega,(BU)) + bar(N)[t]^list(omega,(BL)))
      )
      },
    "s" = {
      METHOD <- "Foster-Stuart s-statistic"
      ylabel <- ifelse(all(w == 1), 
                       "Mean s-statistic",
                       "Mean weighted s-statistic"
                       #expression(bar(N)[t]^list((FU)) + bar(N)[t]^list((FL))),
                       #expression(bar(N)[t]^list(omega,(FU)) + bar(N)[t]^list(omega,(FL)))
      )
      },
    "U" = {
      METHOD <- "Forward - backward upper records"
      ylabel <- ifelse(all(w == 1), 
                       "Mean U-statistic",
                       "Mean weighted U-statistic"
                       #expression(bar(N)[t]^list((FU)) - bar(N)[t]^list((BU))),
                       #expression(bar(N)[t]^list(omega,(FU)) - bar(N)[t]^list(omega,(BU)))
      )
      },
    "L" = {
      METHOD <- "Backward - forward lower records"
      ylabel <- ifelse(all(w == 1), 
                       "Mean L-statistic",
                       "Mean weighted L-statistic"
                       #expression(bar(N)[t]^list((BL)) - bar(N)[t]^list((FL))),
                       #expression(bar(N)[t]^list(omega,(BL)) - bar(N)[t]^list(omega,(FL)))
      )
      },
    "W" = {
      METHOD <- "Diersen-Trenkler W-statistic"
      ylabel <- ifelse(all(w == 1), 
                       "Mean W-statistic",
                       "Mean weighted W-statistic"
                       #expression(bar(N)[t]^list((FU)) + bar(N)[t]^list((BL))),
                       #expression(bar(N)[t]^list(omega,(FU)) + bar(N)[t]^list(omega,(BL)))
      )
      })
  fun   <- deparse(weights)[2]
  
  if (!all(w == 1)) { METHOD <- paste(METHOD, "with weights =", fun) }
  ###################################  
  df <- data.frame(t)
  
  n <- ifelse(statistic %in% c("s", "I"), w[1] * 2, 0) 
  df[2:3] <- matrix(c(n, n, sapply(2:Trows, FUN = function(iter) .E.foster.plot(X[1:iter,], weights, statistic))), ncol = 2, byrow = TRUE)
  colnames(df)[2:3] <- c("stat", "mu")
  
  if (conf.int) {
    stat <- mu <- CI1 <- CI2 <- NULL
    df$sigma <- c(0, sqrt(sapply(2:Trows, FUN = .VAR.foster.plot, weights, statistic) / Mcols))
    q <- stats::qnorm((1 - conf.level) / 2)
    df$CI1 <- df$mu + q * df$sigma   #lower bound
    df$CI2 <- df$mu - q * df$sigma   #upper bound
  }
  
  # plot
  graf <- 
    ggplot2::ggplot(data = df, ggplot2::aes(x = t, y = stat)) + 
    ggplot2::theme_bw() + 
    ggplot2::labs(title = METHOD, subtitle = paste('Data:', DNAME), x = "Time", y = ylabel) +
    ggplot2::geom_line(ggplot2::aes(x = t, y = mu, linetype = "CI"), colour = conf.col) +
    ggplot2::theme(legend.position = 'bottom') +
    ggplot2::scale_color_manual(name = "Statistic", values = c("points" = point.col), label = c("points" = statistic)) +
    ggplot2::scale_shape_manual(name = "Statistic", values = c("points" = point.shape), label = c("points" = statistic)) +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 1), shape = ggplot2::guide_legend(order = 1),linetype = ggplot2::guide_legend(order = 2))
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
  graf <- graf +
    ggplot2::geom_point(ggplot2::aes(colour = "points", shape = "points"))
  
  return(graf)
}

.E.foster.plot <- function(X, 
                        weights,
                        statistic) {

  Trows <- NROW(X)
  t <- 1:Trows
  w <- weights(t)

  switch(statistic,
         "D" = {
           Xrev <- series_rev(X)
           I <- .I.record(X, record = "upper", Trows = Trows) - 
             .I.record(   X, record = "lower", Trows = Trows) - 
             .I.record(Xrev, record = "upper", Trows = Trows) + 
             .I.record(Xrev, record = "lower", Trows = Trows) 
           stat <- sum(w * rowMeans(I))
           mu   <- 0
         },
         "d" = {
           I <- .I.record(X, record = "upper", Trows = Trows) -
             .I.record(X, record = "lower", Trows = Trows)
           stat <- sum(w * rowMeans(I))
           mu   <- 0
         },
         "S" = {
           Xrev <- series_rev(X)
           I <- .I.record(X, record = "upper", Trows = Trows) +
             .I.record(   X, record = "lower", Trows = Trows) - 
             .I.record(Xrev, record = "upper", Trows = Trows) - 
             .I.record(Xrev, record = "lower", Trows = Trows) 
           stat <- sum(w * rowMeans(I))
           mu   <- 0
         },
         "s" = {
           I <- .I.record(X, record = "upper", Trows = Trows) + 
            .I.record(X, record = "lower", Trows = Trows)
           stat <- sum(w * rowMeans(I))
           mu   <- 2 * sum(w / t)
         },
         "U" = {
           Xrev <- series_rev(X)
           I <- .I.record(X, record = "upper", Trows = Trows) - 
             .I.record(Xrev, record = "upper", Trows = Trows) 
           stat <- sum(w * rowMeans(I))
           mu   <- 0
         },
         "L" = {
           Xrev <- series_rev(X)
           I <- .I.record(Xrev, record = "lower", Trows = Trows) - 
             .I.record(X, record = "lower", Trows = Trows) 
           stat <- sum(w * rowMeans(I))
           mu   <- 0
         },
         "W" = {
           Xrev <- series_rev(X)
           I <- .I.record(X, record = "upper", Trows = Trows) + 
             .I.record(Xrev, record = "lower", Trows = Trows)
           stat <- sum(w * rowMeans(I))
           mu   <- 2 * sum(w / t)
         }
  )
  
  return(c(stat, mu))
}

.VAR.foster.plot <- function(iter, 
                             weights,
                             statistic) {
  
  Trows <- iter
  t <- 1:Trows
  w <- weights(t)
  
  sigma2 <- switch(statistic,
         "D" = {
           ifelse(length(w) == 1,
                  w^2 * 4 * sum(((t[-1] + 1 + (t[-1] - 1) / choose(Trows, t[-1])) / t[-1]^2)),
                  4 * (.var_N(t, w) - .cov_NFU_NFL(t, w) - .cov_NFU_NBU(t, w) + .cov_NFU_NBL(t, w))
           )
         },
         "d" = { 2 * (.var_N(t, w) - .cov_NFU_NFL(t, w)) },
         "S" = { 4 * (.var_N(t, w) + .cov_NFU_NFL(t, w) - .cov_NFU_NBU(t, w) - .cov_NFU_NBL(t, w)) },
         "s" = { 2 * (.var_N(t, w) + .cov_NFU_NFL(t, w)) },
         "U" = { 2 * (.var_N(t, w) - .cov_NFU_NBU(t, w)) },
         "L" = { 2 * (.var_N(t, w) - .cov_NFU_NBU(t, w)) },
         "W" = { 2 * (.var_N(t, w) + .cov_NFU_NBL(t, w)) }
  )
  
  return(sigma2)
}
