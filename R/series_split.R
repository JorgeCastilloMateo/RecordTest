#' @title Splitted series
#' @description The vector \code{X} of length \eqn{T} is broken into 
#'   \code{Mcols} blocks, each part containing \eqn{T/\code{Mcols}} elements.
#' 
#'   If the vector \code{X} represents consecutive daily values, e.g., and then 
#'   \code{Mcols = 365} is preferred; this function rearranges \code{X} into a 
#'   matrix format, where each column is the vector of values at the same day 
#'   of the year. If it were monthly data it would be preferred, e.g., 
#'   \code{Mcols = 12}.
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the record inference tools in this package.
#'
#'   This function transforms a vector into a matrix, applying the 
#'   following procedure: the first row of the matrix is made up of the first 
#'   \code{Mcols} elements of the vector, the second row by the \code{Mcols} 
#'   following elements, and so on.
#'   The length of the vector must be a multiple of \code{Mcols} 
#'   (see Note otherwise).
#'
#'   In the case of a vector of daily values, \code{Mcols} is usually 365,
#'   so that the first column corresponds to all the values observed at the  
#'   1st of January, the second to the 2nd of January, etc.
#'
#'   If \eqn{X_{i,j}} represents the value in day \eqn{j} of year \eqn{i}, then if
#'   \deqn{\code{X} = (X_{1,1},X_{1,2},\ldots,X_{1,365},X_{2,1},X_{2,2},\ldots,X_{T,365}),}
#'   applying \code{splitSeries} to \code{X} returns the following matrix:
#'   \deqn{
#'   \left( \begin{array}{cccc}
#'        X_{1,1} & X_{1,2} & \cdots & X_{1,365} \\ 
#'        X_{2,1} & X_{2,2} & \cdots & X_{2,365} \\
#'        \vdots & \vdots &  & \vdots \\
#'        X_{T,1} & X_{T,2} & \cdots & X_{T,365}
#'        \end{array} \right)_{T \times 365}. 
#'   }
#'
#' @note \code{\link{series_double}} can be implemented for the same purpose as 
#'   this function but without requiring that the length of \code{X} be 
#'   divisible by \code{Mcols}. It removes the first elements of \code{X} until
#'   its length is divisible by \code{Mcols}.
#' @param X A numeric vector.
#' @param Mcols An integer number, giving the number of columns in the 
#'   final matrix.
#' @return A matrix with \code{Mcols} columns.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_rev}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' series_split(1:100, Mcols = 10)
#' 
#' series_split(TX_Zaragoza$TX)
#'
#' @export series_split
series_split <- function(X, Mcols = 365) {
  
  Trows <- length(X) / Mcols
  
  if(Trows %% 1 != 0) stop("Mcols must be a divisor of the length of X")
  
  return(t(matrix(X, nrow = Mcols, ncol = Trows)))
}