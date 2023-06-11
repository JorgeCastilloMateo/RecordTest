#' @title Record Values
#' 
#' @description Returns the record values of the values in a vector.
#'   A record value is the magnitude of a record observation. 
#'   
#'   If the argument \code{X} is a matrix, then each column is treated as a 
#'   different vector.
#'   
#' @details The sequence of record values \eqn{\{R_1,\ldots,R_I\}} can be 
#'   expressed in terms of the record times 
#'   \code{\link{L.record}} by
#'   \deqn{R_i = X_{L_i}.}
#'   
#'   Record values can be calculated for both upper and lower records.
#'
#' @inheritParams I.record
#' @return If \code{X} is a vector, the function returns a list containing the 
#'   vector of record values. If \code{X} is a matrix, the function returns a 
#'   list where each element is a vector indicating the record values of the 
#'   corresponding \code{X} column.
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{L.record}}, 
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}, 
#'   \code{\link{p.record}}, \code{\link{R.record}}, 
#'   \code{\link{records}}, \code{\link{S.record}}
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
#' R.record(Y1)
#' R.record(Y)
#' 
#' @export R.record
R.record <- function(X, record = c("upper", "lower"), weak = FALSE){
  
  X <- as.matrix(X)
  L <- L.record(X, record = record, weak = weak)
  
  return(apply(as.matrix(1:ncol(X)), 1, function(i) X[L[[i]], i], simplify = FALSE))
}
