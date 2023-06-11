#' @title Record Times
#' 
#' @description Returns the record times of the values in a vector.
#'   The record times are the positions in a vector where a record occurs. 
#'   
#'   If the argument \code{X} is a matrix, then each column is treated as a 
#'   different vector.
#'   
#' @details The sequence of record times \eqn{\{L_1,\ldots,L_I\}} can be 
#'   expressed in terms of the record indicator random variables 
#'   \code{\link{I.record}} by
#'   \deqn{L_i = \min\{ t \mid I_1 + I_2 + \ldots + I_t = i \}.}
#'   
#'   Record times can be calculated for both upper and lower records.
#'
#' @inheritParams I.record
#' @return If \code{X} is a vector, the function returns a list containing the 
#'   vector of record times. If \code{X} is a matrix, the function returns a 
#'   list where each element is a vector indicating the record times of the 
#'   corresponding \code{X} column.
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{N.record}},
#'   \code{\link{Nmean.record}}, \code{\link{p.record}}, 
#'   \code{\link{R.record}}, \code{\link{records}}, \code{\link{S.record}}
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. 
#' Wiley Series in Probability and Statistics. Wiley, New York.
#' \doi{10.1002/9781118150412}.
#' 
#' @examples
#' Y1 <- c( 1,  5,  3,  6,  6,  9,  2)
#' Y2 <- c(10,  5,  3,  6,  6,  9,  2)
#' Y3 <- c( 5,  7,  3,  6, 19,  2, 20)
#' Y  <- cbind(Y1, Y2, Y3)
#' 
#' L.record(Y1)
#' L.record(Y)
#' 
#' @export L.record
L.record <- function(X, record = c("upper", "lower"), weak = FALSE){
  
  X <- I.record(X, record = record, weak = weak)
  
  return(apply(X, 2, function(x) which(x == 1), simplify = FALSE))
}
