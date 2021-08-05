#' @title From Record Times to Time Series
#' @description If only the record times are available (upper or lower, or 
#'   both) and not the complete series, \code{series_record} builds a complete 
#'   series with the same record occurrence as specified in the arguments. 
#'   This function is useful to apply the plots and tests within
#'   \code{\link{RecordTest-package}} to a vector of record times.
#' @param L_upper,L_lower A vector of (increasing) integers denoting the upper 
#'   or/and lower record times.
#' @param R_upper,R_lower (Optional) A vector of (increasing/decreasing)
#'   values denoting the upper or/and lower record values.
#' @param Trows Integer indicating the actual length of the series. If it is 
#'   not specified, then the length of the series is assumed equal to the last
#'   record occurrence.
#' @return A vector of length \code{Trows} with \code{L_upper} upper or/and 
#'   \code{L_lower} lower record times and \code{R_upper} upper or/and 
#'   \code{R_lower} lower record values.
#' @note Remember that the first observation in a series is always a record 
#'   time.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_rev}}, 
#'   \code{\link{series_split}}, \code{\link{series_ties}},
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' # upper record times observed in a 100 length time series
#' L <- c(1, 4, 14, 40, 45, 90)
#' X <- series_record(L_upper = L, Trows = 100)
#' 
#' # now you can apply plots and tests for upper records to the X series
#' #N.plot(X)
#' #N_normal.test(X)
#' 
#' # if you also have lower record times
#' L_lower <- c(1, 2, 12, 56, 57, 78, 91)
#' X <- series_record(L_upper = L, L_lower = L_lower, Trows = 100)
#' 
#' # now you can apply plots and tests to the X series with both types of record times
#' #foster.plot(X, statistic = 'd')
#' #foster.test(X, statistic = 'd')
#' 
#' # apply to the 200-meter Olympic records from 1900 to 2020
#' or200m <- series_record(L_lower = Olympic_records_200m$t, 
#'                         R_lower = Olympic_records_200m$value,
#'                         Trows = 27)
#' # some plots and tests                    
#' N.plot(or200m, record = c(0,1,0,0))                         
#' N.test(or200m, record = "lower", distribution = "poisson-binomial")
#' @export series_record
series_record <- function(L_upper, R_upper, L_lower, R_lower, Trows = NA) {
  
  mu <- missing(L_upper)
  ml <- missing(L_lower)
  muR <- missing(R_upper)
  mlR <- missing(R_lower)
  
  if (!mu) {
    if (!all(L_upper %% 1 == 0)) { stop("'L_upper' should be a vector of integers") }
    if (!all(diff(L_upper) > 0)) { stop("'L_upper' should be an strictly increasing vector") }
    if (L_upper[1] != 1) { stop("1 is not in 'L_upper' and 1 is always a record time") }
    if (!is.na(Trows) && max(L_upper) > Trows) { stop("'Trows' should be bigger than the last record time") }
    if (!muR) {
      if (length(L_upper) != length(R_upper)) { stop("'L_upper' and 'R_upper' should have the same length") }
      if (!all(diff(R_upper) > 0)) { stop("'R_upper' should be an strictly increasing vector") }
    } 
  }
  
  if (!ml) {
    if (!all(L_lower %% 1 == 0)) { stop("'L_lower' is not a vector of integers") }
    if (!all(diff(L_lower) > 0)) { stop("'L_lower' is not an strictly increasing vector") }
    if (L_lower[1] != 1) { stop("1 is not in 'L_lower' and 1 is always a record time") }
    if (!is.na(Trows) && max(L_lower) > Trows) { stop("'Trows' must be bigger than the last record time") }
    if (!mlR) {
      if (length(L_lower) != length(R_lower)) { stop("'L_lower' and 'R_lower' should have the same length") }
      if (!all(diff(R_lower) < 0)) { stop("'R_lower' should be an strictly decreasing vector") }
    }
  }
  
  if (!mu && !ml) {
    
    if (length(intersect(L_upper, L_lower)) != 1) { stop("The only common element between 'L_upper' and 'L_lower' should be 1") }
    if (!muR && !mlR && muR[1] != mlR[1]) { stop("First record value should be the same for 'R_upper' and 'R_lower'") }
    
    if (is.na(Trows)) Trows <- max(L_upper, L_lower)
    
    I <- c(0, diff(series_record(L_upper = L_upper, Trows = Trows)) + 
      diff(series_record(L_lower = L_lower, Trows = Trows)))
    
    if (!muR && !mlR) {
      X <- rep(R_upper[1], Trows)
      X[I ==  1] <- R_upper[-1]
      X[I == -1] <- R_lower[-1]
    } else {
      X <- rep(0, Trows)
      X[I ==  1] <- 1:(length(L_upper) - 1)
      X[I == -1] <- -(1:(length(L_lower) - 1))
    }
  } else {
    
    if (mu && ml) { stop("Argument 'L_upper' or 'L_lower' is missing") }
    
    if (mu) L <- L_lower
    else    L <- L_upper 

    N <- length(L)
    
    if (is.na(Trows)) Trows <- L[N]
    
    X <- rep(1:N, times = c(diff(L), Trows - L[N] + 1)) 
    
    if (!mu && !muR)      X <- R_upper[X]
    else if (!ml && !mlR) X <- R_lower[X]
    else if (mu) X <- -X
  }
    
  return(X)                
}
