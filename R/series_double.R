#' @title Double the Number of Series
#' @description This function changes the format of a matrix transforming a
#'   \eqn{T \times M} matrix in a 
#'   \eqn{\lfloor T/\code{k} \rfloor \times \code{k}\,M} matrix in the 
#'   following way.
#'   
#'   First, the matrix is divided into \code{k} matrices 
#'   \eqn{\lfloor T/\code{k} \rfloor \times M}, containing the rows whose 
#'   remainder of the division of the row number by \code{k} is 
#'   \eqn{1,2,\ldots,\code{k}-1,0}, respectively;
#'   and secondly those matrices are \code{cbind}ed.
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the exploratory and inference tools based on 
#'   theory of records within this package.
#'
#'   Most of the record inference tools require a high number of independent 
#'   series \eqn{M} (number of columns) to be applied.
#'   If \eqn{M} is low and the time period of observation, \eqn{T}, is high 
#'   enough, the following procedure can be applied in order to multiply by
#'   \code{k} the value \eqn{M}.
#'   The approach  consists of considering that the observations at two 
#'   (or more) consecutive times, \eqn{t} and \eqn{t+1} (or \eqn{t+\code{k}-1}),
#'   are independent observations measured at the same time unit.
#'   That means that we are doubling (or multiplying by \code{k}) the original 
#'   time unit of the records, so that the length of the observation period 
#'   will be \eqn{\lfloor T/\code{k} \rfloor}.
#'   This function rearranges the original data matrix into the new format.
#'
#'   If the number of rows of the original matrix is not divisible by \code{k},
#'   the first \code{nrow(X) \%\% k} rows are deleted.
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param k Integer \eqn{> 1}, times to increase the number of columns.
#' @return A \eqn{\lfloor T/\code{k} \rfloor \times \code{k}\,M} matrix.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_record}}, \code{\link{series_rev}},
#'   \code{\link{series_split}}, \code{\link{series_ties}}
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' series_double(matrix(1:100, 10, 10))
#' 
#' series_double(ZaragozaSeries, k = 4)
#' 
#' @export series_double
series_double <- function(X, k = 2) {
  
  X <- as.matrix(X)
  Trows <- nrow(X)
  
  if (Trows < k) { stop(paste("'X' does not have enough elements for 'k' =", k)) }
  
  remainder <- Trows %% k
  if (remainder != 0) X <- X[(remainder + 1):Trows, , drop = FALSE]
    
  return(series_split(c(t(X)), Mcols = ncol(X) * k))
}