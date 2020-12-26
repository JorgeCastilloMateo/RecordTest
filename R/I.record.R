#' @title Record Indicators
#' @description Returns the sample record indicators of the values in a vector.
#'   The record indicator for each value in a vector is a binary variable which
#'   takes the value 1 if the corresponding value in the vector is a record and
#'   0 otherwise. 
#'   
#'   If the argument \code{X} is a matrix, then each column is treated as a 
#'   different vector.
#' @details Let \eqn{\{X_1,\ldots,X_T\}} be a vector of random variables of 
#'   size \eqn{T}. An observation \eqn{X_t} will be called an upper record 
#'   value if its value exceeds that of all previous observations. An 
#'   analogous definition deals with lower record values.
#'   Here, \eqn{X_1} is referred to as the reference value or the trivial record.
#'   Then, the sequence of record indicator random variables 
#'   \eqn{\{I_1,\ldots,I_T\}} is given by
#'       \deqn{I_t = \left\{ 
#'         \begin{array}{ll} 
#'           1 & \mbox{if } X_t \mbox{ is a record,} \\ 
#'           0 & \mbox{if } X_t \mbox{ is not a record.} 
#'         \end{array} \right.} 
#'
#'   The method \code{I.record} calculates the sample sequence above if the 
#'   argument \code{X} is a numeric vector. If the argument \code{X} is a 
#'   matrix (or data frame) with \eqn{M} columns, the method \code{I.record} 
#'   calculates the sample sequence above for each column of the object as if 
#'   all columns were different sequences.
#'  
#'   Summarily:
#'   \deqn{\code{I.record}: \code{X} = \left(
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
#'   All the procedure above can be extended to weak records, which also count
#'   the ties as a new (weak) record. Ties are possible in discrete variables
#'   or if a continuous variable has been rounded. Weak records can be computed
#'   if \code{weak = TRUE}.
#'  
#' @aliases I.record.default I.record.numeric I.record.matrix I.record
#' @param X A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be 
#'   calculated, "upper" or "lower".
#' @param weak Logical. If \code{TRUE}, weak records are also counted. Default
#'   to \code{FALSE}.
#' @return A binary matrix of the same length or dimension as \code{X}, 
#'   indicating the record occurrence.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{L.record}}, 
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}, 
#'   \code{\link{p.record}}, \code{\link{R.record}},
#'   \code{\link{records}}, \code{\link{S.record}}
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. 
#' Wiley Series in Probability and Statistics. Wiley, New York.
#' @examples
#' X <- c(1, 5, 3, 6, 6, 9, 2, 11, 17, 8)
#' I.record(X)
#' I.record(X, weak = TRUE)
#' 
#' I.record(ZaragozaSeries)
#' # record argument can be shortened
#' I.record(ZaragozaSeries, record = "l")
#' 
#' @export I.record
I.record <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  UseMethod("I.record", X)
}

#' @rdname I.record
#' @method I.record default
#' @export 
I.record.default <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  return(I.record.matrix(X = as.matrix(X), record = record, weak = weak))
}

#' @rdname I.record
#' @method I.record numeric
#' @export 
I.record.numeric <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  record       <- match.arg(record)
  record_upper <- record == "upper"
  
  if (record_upper & !weak) {
    I <- c(1, cummax(X)[-length(X)] < X[-1]) 
  } else if (!record_upper & !weak) {
    I <- c(1, cummin(X)[-length(X)] > X[-1])
  } else if (record_upper & weak) {
    I <- c(1, cummax(X)[-length(X)] <= X[-1]) 
  } else { # !record_upper & weak
    I <- c(1, cummin(X)[-length(X)] >= X[-1]) 
  } 
  
  names(I) <- names(X) 
    
  return(as.matrix(I))
}

#' @rdname I.record
#' @method I.record matrix
#' @export 
I.record.matrix <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  record <- match.arg(record)
  
  DNAME <- rownames(X)
  X <- apply(X, 2, I.record.numeric, record = record, weak = weak)
  rownames(X) <- DNAME
  
  return(X)
}

##########################
### INTERNAL FUNCTIONS ###
##########################
.I.record <- function(X, ...) { UseMethod(".I.record") }

#' @method .I.record default
.I.record.default <- function(X, record, Trows) {
  
  return(.I.record.matrix(X = as.matrix(X), record = record, Trows = Trows))
}

#' @method .I.record numeric
.I.record.numeric <- function(X, record, Trows) {
  
  if (record == "upper") { I <- c(1, cummax(X)[-length(X)] < X[-1]) }
  else                   { I <- c(1, cummin(X)[-length(X)] > X[-1]) }
  
  return(as.matrix(I))
}

#' @method .I.record matrix
.I.record.matrix <- function(X, record, Trows) {
  
  if (record == "upper") { I <- apply(X, 2, .I.upper, Trows = Trows) }
  else                   { I <- apply(X, 2, .I.lower, Trows = Trows) }
  
  return(I)
}

.I.upper <- function(X, Trows) { return(c(1, cummax(X)[-Trows] < X[-1])) }
.I.lower <- function(X, Trows) { return(c(1, cummin(X)[-Trows] > X[-1])) }
