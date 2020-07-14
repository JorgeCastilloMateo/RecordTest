#' @title Number and Mean Number of Records up to time \eqn{t}
#' @description \code{N.record} and \code{Nmean.record} calculates
#'   the number of records up to time \eqn{t} in a vector (or every matrix 
#'   column) of length \eqn{T} and the mean number of records up to time \eqn{t}
#'   in \eqn{M} vectors (or matrix columns), respectively.
#' @details The record counting process \eqn{\{N_1,\ldots,N_T\}} is defined by 
#'   the number of records up to time \eqn{t}, and can be expressed in terms of 
#'   the record indicator random variables \code{\link{I.record}} by
#'   \deqn{N_t = I_1 + I_2 + \ldots + I_t.}
#'
#'   If \code{XM_T} is a matrix with \eqn{M>1} columns, each column is treated 
#'   as a vector and \code{Nmean.record} calculates for each \eqn{t},
#'   \deqn{\bar N_t = \frac{N_{t,1}+ \ldots + N_{t,M}}{M}.}
#'   
#'   Summarily:
#'   \deqn{\code{N.record}: \code{XM\_T} = \left(
#'                  \begin{array}{cccc} 
#'                    X_{1,1} & X_{1,2} & \cdots & X_{1,M} \\ 
#'                    X_{2,1} & X_{2,2} & \cdots & X_{2,M} \\ 
#'                    \vdots & \vdots &  & \vdots \\ 
#'                    X_{T,1} & X_{T,2} & \cdots & X_{T,M} \\ 
#'                  \end{array} \right) 
#'                  \longrightarrow
#'                  \left(
#'                  \begin{array}{cccc} 
#'                    N_{1,1} & N_{1,2} & \cdots & N_{1,M} \\ 
#'                    N_{2,1} & N_{2,2} & \cdots & N_{2,M} \\ 
#'                    \vdots & \vdots &  & \vdots \\ 
#'                    N_{T,1} & N_{T,2} & \cdots & N_{T,M} \\ 
#'                  \end{array} \right)}
#'   and
#'   \deqn{\code{Nmean.record}: \code{XM\_T} 
#'         \longrightarrow
#'         \Big( \bar{N}_1, \bar{N}_2, \cdots, \bar{N}_T \Big).}
#'         
#'   Number and mean number of records for both upper and lower records can be 
#'   calculated.
#'   
#' @note If \code{XM_T} is a vector both functions return the same values, 
#'   \code{N.record} as a matrix and \code{Nmean.record} as a vector.
#' @aliases N.record Nmean.record
#' @inheritParams I.record
#' @return \code{N.record} returns a numeric matrix with the number of records 
#'   up to each time (row) \eqn{t} for a vector or each column in \code{XM_T}. 
#'   \code{Nmean.record} returns a numeric vector with the mean number of 
#'   records up to each time (row) \eqn{t}. 
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{L.record}}, 
#'   \code{\link{M.record}}, \code{\link{P.record}}, \code{\link{records}}
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). \emph{Records}. 
#' New York: Wiley.
#' @examples
#' Y1 <- c( 1,  5,  3,  6,  6,  9,  2)
#' Y2 <- c(10,  5,  3,  6,  6,  9,  2)
#' Y3 <- c( 5,  7,  3,  6, 19,  2, 20)
#' Y  <- cbind(Y1, Y2, Y3)
#' 
#' N.record(Y)
#' Nmean.record(Y)
#' 
#' N.record(ZaragozaSeries)
#' Nmean.record(ZaragozaSeries, record = 'l')
#' 
#' @export N.record
N.record <- function(XM_T, record = c('upper', 'lower')) {
  
  XM_T <- I.record(XM_T, record = record)
  
  return(apply(XM_T, 2, cumsum))
}

#' @rdname N.record
#' @export Nmean.record
Nmean.record <- function(XM_T, record = c('upper', 'lower')) {
  
  XM_T <- N.record(XM_T, record = record)
  
  return(rowMeans(XM_T))
}

