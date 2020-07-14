#' @title Plot of the record times
#' @importFrom ggplot2 ggplot aes theme_bw theme element_blank element_rect
#'   element_text labs geom_hline geom_point xlim facet_grid vars 
#' @description This function constructs a ggplot object to display the upper 
#'   and lower record times for both forward and backward sequences.
#' @details The function can be applied to plot the record times in a vector 
#'   (if argument \code{XM_T} is a vector) or to plot and compare the record 
#'   times in a set of vectors (if argument \code{XM_T} is a matrix). In the 
#'   latter case, the approach to obtain the record times is applied to each 
#'   column of the matrix.
#'   
#'   A matrix of four panels is displayed for upper and lower records, and for
#'   the forward and backward (\code{\link{series_rev}}) vectors.
#'
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param colour_point Colour to plot points.
#' @param colour_line Colour to plot lines.
#' @return A ggplot object.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L.record}}, \code{\link{L_lr.test}}
#' @examples
#' Y <- c(1, 5, 3, 6, 6, 9, 2, 11, 17, 8)
#' L.plot(Y)
#' 
#' L.plot(ZaragozaSeries, colour_point = 1)
#'
#' @export L.plot

L.plot <- function(XM_T, colour_point = 'skyblue3', colour_line = 'grey95'){
  
  DNAME <- deparse(substitute(XM_T))
  XM_T  <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  Mcols <- ncol(XM_T)
  
  XM_Trev <- apply(XM_T, 2, rev)
  
  I     <- I.record(XM_T, record = 'upper')
  IL    <- I.record(XM_T, record = 'lower')
  Irev  <- I.record(XM_Trev, record = 'upper')
  ILrev <- I.record(XM_Trev, record = 'lower')
  
  L.fun <- function(I2) {
    
    L <- cbind(which(I2[,1] == 1), 1)
  
    n <- 1
    while (n < Mcols) {
      n <- n + 1
      L <- rbind(L, cbind(which(I2[,n] == 1), n))
    }
    return(L)
  }
  
  L     <- L.fun(I)
  LL    <- L.fun(IL)
  Lrev  <- L.fun(Irev)
  LLrev <- L.fun(ILrev)
  
  L     <- data.frame(L, roundTrip = rep('forward', nrow(L)), upperLower = rep('upper', nrow(L)))
  LL    <- data.frame(LL, roundTrip = rep('forward', nrow(LL)), upperLower = rep('lower', nrow(LL)))
  Lrev  <- data.frame(Lrev, roundTrip = rep('backward', nrow(Lrev)), upperLower = rep('upper', nrow(Lrev)))
  LLrev <- data.frame(LLrev, roundTrip = rep('backward', nrow(LLrev)), upperLower = rep('lower', nrow(LLrev)))
  
  L <- rbind(L, LL, Lrev, LLrev)
  
  # plot
  graf <- 
    ggplot2::ggplot(data = L, ggplot2::aes(x = L[,1], y = L[,2])) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          strip.background = ggplot2::element_rect(fill = 'gray23'),
          strip.text = ggplot2::element_text(colour = 'white')) +
    ggplot2::labs(title = "Record times", 
         subtitle = paste('Data:', DNAME),
         x = "Times", y = "Series") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = L[,2]), colour = colour_line) + 
    ggplot2::geom_point(colour = colour_point) + ggplot2::xlim(0, Trows) +
    ggplot2::facet_grid(ggplot2::vars(L$roundTrip), ggplot2::vars(L$upperLower))
  ###################################
  
  return(graf)
}