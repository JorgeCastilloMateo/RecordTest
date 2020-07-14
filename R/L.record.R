#' @title Record Times
#' @description This function calculates the times (positions in the vector) 
#'   where records occur.
#' @details If \code{XM_T} is a matrix, the approach to obtain record times is 
#'   applied to each column of the matrix. 
#'   
#'   Record times can be calculated for both upper and lower records.
#'
#' @inheritParams I.record
#' @return If \code{XM_T} is a vector, the function returns a column matrix 
#'   containing the record times. If \code{XM_T} is a matrix, the function 
#'   returns a list where each element is a vector indicating the record times 
#'   of the corresponding \code{XM_T} column.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{M.record}},
#'   \code{\link{N.record}}, \code{\link{P.record}}, \code{\link{Nmean.record}},
#'   \code{\link{records}}
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
L.record <- function(XM_T, record = c('upper', 'lower')){
  
  XM_T <- I.record(XM_T, record = record)
  
  return(apply(XM_T, 2, function(x) which(x == 1)))
}
