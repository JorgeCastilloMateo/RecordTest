#' @title Times of Record Plot
#' @importFrom ggplot2 ggplot aes theme_bw theme element_blank element_rect
#'   element_text labs geom_hline geom_point xlim facet_grid vars 
#' @description This function builds a ggplot object to display the upper 
#'   and lower record times for both forward and backward directions.
#' @details The function can be applied to plot the record times in a vector 
#'   (if argument \code{X} is a vector) or to plot and compare the record 
#'   times in a set of vectors (if argument \code{X} is a matrix). In the 
#'   latter case, the approach to obtain the record times is applied to each 
#'   column of the matrix.
#'   
#'   If \code{all = TRUE}, a matrix of four panels is displayed for upper and
#'   lower records, and for the forward and backward (\code{\link{series_rev}})
#'   directions. Otherwise, only one type of forward record is displayed.
#'   
#'   An example of use of a plot with similar ideas is shown in Benestad 
#'   (2004, Figures 3 and 8).
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param all Logical. If \code{TRUE} (the default) the four types of record
#'   are displayed.
#' @param record If \code{all = FALSE}, a character string indicating the type
#'   of record to be calculated, "upper" or "lower".
#' @param point.col,point.alpha Color and transparency of the points.
#' @param line.col Color to plot lines.
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L.record}}
#' @references
#' Benestad RE (2004). “Record-Values, Nonstationarity Tests and Extreme Value Distributions.”
#' \emph{Global and Planetary Change}, \strong{44}(1-4), 11-26.
#' \href{https://doi.org/10.1016/j.gloplacha.2004.06.002}{doi:10.1016/j.gloplacha.2004.06.002}
#' 
#' @examples
#' Y <- c(1, 5, 3, 6, 6, 9, 2, 11, 17, 8)
#' L.plot(Y, all = FALSE)
#' 
#' L.plot(ZaragozaSeries, point.col = 1)
#'
#' @export L.plot

L.plot <- function(X, 
                   all = TRUE, 
                   record = c("upper", "lower"), 
                   point.col = "gray23", 
                   point.alpha = 0.8, 
                   line.col = "gray95"){
  
  record <- match.arg(record)
  
  DNAME <- deparse(substitute(X))
  
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  if (Trows == 1) { stop("'NROW(X)' should be greater than 1") }

  L.fun <- function(I) {
    
    L <- apply(I, 2, FUN = function(x) which(x == 1))
    
    if (!is.list(L)) L <- list(L)
    
    return(cbind(unlist(L), rep(seq_len(Mcols), lengths(L))))
  }
  
  if (all) {
    
    Xrev <- series_rev(X)
    
    I     <- .I.record(X,    record = "upper", Trows = Trows)
    IL    <- .I.record(X,    record = "lower", Trows = Trows)
    Irev  <- .I.record(Xrev, record = "upper", Trows = Trows)
    ILrev <- .I.record(Xrev, record = "lower", Trows = Trows)
    
    L     <- L.fun(I)
    LL    <- L.fun(IL)
    Lrev  <- L.fun(Irev)
    LLrev <- L.fun(ILrev)
    
    L     <- data.frame(L, roundTrip = rep("forward", nrow(L)), upperLower = rep("upper", nrow(L)))
    LL    <- data.frame(LL, roundTrip = rep("forward", nrow(LL)), upperLower = rep("lower", nrow(LL)))
    Lrev  <- data.frame(Lrev, roundTrip = rep("backward", nrow(Lrev)), upperLower = rep("upper", nrow(Lrev)))
    LLrev <- data.frame(LLrev, roundTrip = rep("backward", nrow(LLrev)), upperLower = rep("lower", nrow(LLrev)))
    
    L <- rbind(L, LL, Lrev, LLrev)
    
    # plot
    graf <- 
      ggplot2::ggplot(data = L, ggplot2::aes(x = L[,1], y = L[,2])) + 
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     strip.background = ggplot2::element_rect(fill = "black"),
                     strip.text = ggplot2::element_text(colour = "white")) +
      ggplot2::labs(title = "Times of record", 
                    subtitle = paste("Data:", DNAME),
                    x = "Time", y = "Series") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = L[,2]), colour = line.col) + 
      ggplot2::geom_point(colour = point.col, alpha = point.alpha) + ggplot2::xlim(0, Trows) +
      ggplot2::facet_grid(ggplot2::vars(L$roundTrip), ggplot2::vars(L$upperLower))
    ###################################
    
  } else {
    
    L <- L.fun(I.record(X, record = record))

    # plot
    graf <- 
      ggplot2::ggplot(data = data.frame(L), ggplot2::aes(x = L[,1], y = L[,2])) + 
      ggplot2::theme_bw() + 
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::labs(title = paste("Times of", ifelse(record == "upper", "upper", "lower"), "record"), 
                    subtitle = paste("Data:", DNAME),
                    x = "Time", y = "Series") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = L[,2]), colour = line.col) + 
      ggplot2::geom_point(colour = point.col, alpha = point.alpha) + ggplot2::xlim(0, Trows)
    ###################################
  }
  
  return(graf)
}