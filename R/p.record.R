#' @title Probabilities of Record
#' @description 
#'   \code{S.record} and \code{p.record} return the sample number of 
#'   records and mean number of records at each time \eqn{t} in a set of \eqn{M}
#'   vectors (columns of \code{X}), respectively. In particular, 
#'   \code{p.record} is the estimated record probability at each time \eqn{t}.
#'   
#'   (For the introduccion to records see Details in \code{\link{I.record}}.) 
#' @details Given a matrix formed by \eqn{M} vectors (columns), measured at 
#'   \eqn{T} times (rows), \code{M.record} calculates the number of records in 
#'   the \eqn{M} vectors at each observed time \eqn{t}, \eqn{S_t}.
#'
#'   The function \code{p.record} is equivalent, but calculates the proportion 
#'   of records at each time \eqn{t}, that is the ratio:
#'   \deqn{\hat p_t = \frac{S_t}{M} = \frac{I_{t,1} + \ldots + I_{t,M}}{M},}
#'   this proportion is an estimation of the probability of record at that time.
#'  
#'   Following the notation in \code{\link{I.record}}, in summary:
#'   \deqn{\code{X} = \left(
#'                  \begin{array}{cccc} 
#'                    X_{1,1} & X_{1,2} & \cdots & X_{1,M} \\ 
#'                    X_{2,1} & X_{2,2} & \cdots & X_{2,M} \\ 
#'                    \vdots & \vdots &  & \vdots \\ 
#'                    X_{T,1} & X_{T,2} & \cdots & X_{T,M} \\ 
#'                  \end{array} \right) 
#'                  \begin{array}{lc} 
#'                  \stackrel{\code{S.record}}{\longrightarrow} &
#'                  \Big( S_1, S_2, \cdots, S_T \Big) \\ \\ 
#'                  \stackrel{\code{p.record}}{\longrightarrow} &
#'                  \Big( \hat p_1, \hat p_2, \cdots, \hat p_T \Big) \\
#'                  \end{array}}
#' 
#'   Summaries for both upper and lower records can be calculated.
#' 
#' @aliases S.record p.record
#' @inheritParams I.record
#' @return A vector with the number (or proportion in the case of 
#'   \code{p.record}) of records at each time \eqn{t} (row).
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{L.record}}, 
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}, 
#'   \code{\link{R.record}}, \code{\link{records}}
#' @references 
#' Cebrián A, Castillo-Mateo J, Asín J (2021).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' Available at Research Square \doi{10.21203/rs.3.rs-214787/v1}
#' 
#' @examples
#' Y1 <- c( 1,  5,  3,  6,  6,  9,  2)
#' Y2 <- c(10,  5,  3,  6,  6,  9,  2)
#' Y3 <- c( 5,  7,  3,  6, 19,  2, 20)
#' Y  <- cbind(Y1, Y2, Y3)
#' 
#' S.record(Y)
#' p.record(Y)
#' 
#' S.record(ZaragozaSeries)
#' p.record(ZaragozaSeries, record = "l")
#' 
#' @export p.record
p.record <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  X <- I.record(X, record = record, weak = weak)
  
  return(rowMeans(X))
}

#' @rdname p.record
#' @export S.record
S.record <- function(X, record = c("upper", "lower"), weak = FALSE) {
  
  X <- I.record(X, record = record, weak = weak)
  
  return(rowSums(X))
}
