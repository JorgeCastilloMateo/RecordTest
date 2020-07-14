#' @title Transforms a \eqn{T \times M} in a 
#'   \eqn{\lfloor T/\code{k} \rfloor \times \code{k}\,M} matrix
#' @description This function changes the format of a matrix, 
#'   in the following way.
#'   First, the matrix is divided into \code{k} matrices 
#'   \eqn{\lfloor T/\code{k} \rfloor \times M}, containing the rows whose 
#'   remainder of the division of the row number by \code{k} is 
#'   \eqn{1,2,\ldots,\code{k}-1,0}, respectively;
#'   and secondly those matrices are cbinded.
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the record inference tools in this package.
#'
#'   Most of the record inference tools require a high number of independent 
#'   series \eqn{M} (number of columns) to be applied.
#'   If \eqn{M} is low and the time period of observation, \eqn{T}, is high 
#'   enough, the following procedure can be applied in order to multiply by
#'   \code{k} the value \eqn{M}.
#'   The approach  consists of considering that the observations at two 
#'   (or more) consesutive times, \eqn{t} and \eqn{t+1} (or \eqn{t+\code{k}-1}),
#'   are independent observations measured at the same time unit.
#'   That means that we are doubling (or multiplying by \code{k}) the original 
#'   time unit of the records, so that the length of the observation period 
#'   will be \eqn{\lfloor T/\code{k} \rfloor}.
#'   This function rearranges the original data matrix into the new format.
#'
#'   If the number of rows of the original matrix is not divisible by \code{k},
#'   the first \code{nrow(XM_T) \%\% k} rows are deleted.
#'
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param k Integer \eqn{> 1}, times to increase the number of columns.
#' @return A \eqn{\lfloor T/\code{k} \rfloor \times \code{k}\,M} matrix.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' series_double(matrix(1:100, 10, 10))
#' 
#' series_double(ZaragozaSeries, k = 4)
#' 
#' @export series_double
series_double <- function(XM_T, k = 2) {
  
  XM_T <- as.matrix(XM_T)
  Trows <- nrow(XM_T)
  
  if (Trows < k) stop(paste("XM_T does not have enough elements for k =", k))
  
  remainder <- Trows %% k
  if (remainder != 0) XM_T <- XM_T[(remainder + 1):Trows, , drop = FALSE]
    
  series_split(c(t(XM_T)), Mcols = ncol(XM_T) * k)
}