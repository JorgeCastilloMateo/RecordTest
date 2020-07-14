#' @title Record Indicator random variables
#' @description This function calculates for each value in the vector a binary 
#'   variable which takes the value 1 if the corresponding value in the vector 
#'   is a record and 0 otherwise. If the argument \code{XM_T} is a matrix, 
#'   each column is treated as a different vector.
#' @details Let \eqn{\{X_1,\ldots,X_T\}} be a sequence of random variables of 
#'   size \eqn{T}. An observation \eqn{X_t} will be called an upper record 
#'   value if its value exceeds that of all previous observations. An 
#'   analogous definition deals with lower record values.
#'   Here \eqn{X_1} is referred to as the reference value or the trivial record.
#'   Then, the sequence of record indicator random variables 
#'   \eqn{\{I_1,\ldots,I_T\}} is given by
#'       \deqn{I_t = \left\{ 
#'         \begin{array}{ll} 
#'           1 & \mbox{if } X_t \mbox{ is a record,} \\ 
#'           0 & \mbox{if } X_t \mbox{ is not a record.} 
#'         \end{array} \right.} 
#'
#'   The method \code{I.record} calculates the sample sequence above if the 
#'   argument \code{XM_T} is a numeric vector. If the argument \code{XM_T} is a 
#'   matrix (or data frame) with \eqn{M} columns, the method \code{I.record} 
#'   calculates the sample sequence above for each column of the object as if 
#'   all columns were different sequences.
#'  
#'   Summarily:
#'   \deqn{\code{I.record}: \code{XM\_T} = \left(
#'                  \begin{array}{cccc} 
#'                    X_{1,1} & X_{1,2} & \cdots & X_{1,M} \\ 
#'                    X_{2,1} & X_{2,2} & \cdots & X_{2,M} \\ 
#'                    \vdots & \vdots &  & \vdots \\ 
#'                    X_{T,1} & X_{T,2} & \cdots & X_{T,M} \\ 
#'                  \end{array} \right) 
#'                  \longrightarrow
#'                  \left(
#'                  \begin{array}{cccc} 
#'                    I_{1,1} & I_{1,2} & \cdots & I_{1,M} \\ 
#'                    I_{2,1} & I_{2,2} & \cdots & I_{2,M} \\ 
#'                    \vdots & \vdots &  & \vdots \\ 
#'                    I_{T,1} & I_{T,2} & \cdots & I_{T,M} \\ 
#'                  \end{array} \right).}
#'                  
#'   Indicators of record occurrence can be calculated for both upper and 
#'   lower records.
#'  
#' @aliases I.record.default I.record.numeric I.record.matrix I.record
#' @param XM_T  A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be 
#'   calculated, "upper" or "lower".
#' @return A binary matrix, indicating the record occurrence.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L.record}}, \code{\link{M.record}}, 
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}, 
#'   \code{\link{P.record}}, \code{\link{records}}
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. New York: Wiley.
#' @examples
#' X <- c(1, 5, 3, 6, 6, 9, 2, 11, 17, 8)
#' I.record(X)
#' 
#' I.record(ZaragozaSeries)
#' # record argument can be shortened
#' I.record(ZaragozaSeries, record = 'l')
#' 
#' @export I.record
I.record <- function(XM_T, record = c('upper', 'lower')) {
  
  UseMethod('I.record', XM_T)
}

#' @rdname I.record
#' @method I.record default
#' @export 
I.record.default <- function(XM_T, record = c('upper', 'lower')) {
  
  I.record.matrix(XM_T = as.matrix(XM_T), record = record)
}

#' @rdname I.record
#' @method I.record numeric
#' @export 
I.record.numeric <- function(XM_T, record = c('upper', 'lower')) {
  
  record <- match.arg(record)
  
  if (record == 'upper') I <- c(1, cummax(XM_T)[-length(XM_T)] < XM_T[-1]) 
  else                   I <- c(1, cummin(XM_T)[-length(XM_T)] > XM_T[-1])
  
  return(as.matrix(I))
}

#' @rdname I.record
#' @method I.record matrix
#' @export 
I.record.matrix <- function(XM_T, record = c('upper', 'lower')) {
  
  record <- match.arg(record)
  
  Trows <- nrow(XM_T)
  
  if (record == 'upper') I <- apply(XM_T, 2, I.upper, Trows = Trows)
  else                   I <- apply(XM_T, 2, I.lower, Trows = Trows)
  
  return(I)
}

I.upper <- function(XM_T, Trows) return(c(1, cummax(XM_T)[-Trows] < XM_T[-1]))

I.lower <- function(XM_T, Trows) return(c(1, cummin(XM_T)[-Trows] > XM_T[-1]))
