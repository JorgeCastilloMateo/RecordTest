#' @title Split Series
#' 
#' @description The vector \code{X} of length \eqn{T} is broken into 
#'   \code{Mcols} blocks, each part containing \eqn{T/\code{Mcols}} elements.
#' 
#'   If the vector \code{X} represents consecutive daily values in a year, then 
#'   \code{Mcols = 365} is preferred. This function rearranges \code{X} into a 
#'   matrix format, where each column is the vector of values at the same day 
#'   of the year. For monthly data in a year, \code{Mcols = 12} should be used.
#'   
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the exploratory and inference tools based on 
#'   theory of records within this package when the time series presents 
#'   seasonality.
#'
#'   This function transforms a vector into a matrix, applying the 
#'   following procedure: the first row of the matrix is built of the first 
#'   \code{Mcols} elements of the vector, the second row by the \code{Mcols} 
#'   following elements, and so on.
#'   The length of the vector must be a multiple of \code{Mcols} 
#'   (see Note otherwise).
#'
#'   In the case of a vector of daily values, \code{Mcols} is usually 365,
#'   so that the first column corresponds to all the values observed at the  
#'   1st of January, the second to the 2nd of January, etc.
#'
#'   If \eqn{X_{t,m}} represents the value in day \eqn{m} of year \eqn{t}, then if
#'   \deqn{\code{X} = (X_{1,1},X_{1,2},\ldots,X_{1,365},X_{2,1},X_{2,2},\ldots,X_{T,365}),}
#'   applying \code{series_split} to \code{X} returns the following matrix:
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
#'   
#' @param X A numeric vector.
#' @param Mcols An integer number, giving the number of columns in the 
#'   final matrix.
#' @return A matrix with \code{Mcols} columns.
#' 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_record}}, 
#'   \code{\link{series_rev}}, \code{\link{series_ties}},
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @references 
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2), 313-330. 
#' \doi{10.1007/s00477-021-02122-w}.
#' 
#' @examples
#' series_split(1:100, Mcols = 10)
#' 
#' TxZ <- series_split(TX_Zaragoza$TX)
#' dim(TxZ)
#'
#' @export series_split
series_split <- function(X, Mcols = 365) {
  
  Trows <- length(X) / Mcols
  
  if(Trows %% 1 != 0) { stop("'Mcols' must be a divisor of the length of 'X'") }
  
  return(matrix(X, nrow = Trows, ncol = Mcols, byrow = TRUE))
}
