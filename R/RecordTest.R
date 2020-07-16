#' @title RecordTest: A package for testing the classical record model
#' @description RecordTest provides an extensive number of graphical tools 
#'   and tests for randomness to detect trend in location in time series based
#'   on the classical record model. 
#'   It also provides some preprocessing tools often needed to adequately 
#'   prepare a dataset to apply those tests.
#' @details This package provides two families of tests to study the 
#'   hypothesis of the classical record model, that is that the records from a
#'   series of values observed at regular time units (a vector) come from an 
#'   i.i.d. series of continous random variables. A particular problem,
#'   common for example in climate problems, where these tools can be useful, 
#'   is to detect the existence of a trend in a series of variables. If we have
#'   sequences of uncorrelated variables with no seasonal component, the 
#'   hypothesis of i.i.d. variables is equivalent to test the hypothesis of no 
#'   trend.
#'
#'   A sample of \eqn{M} vectors uncorrelated between them is needed to apply 
#'   the implemented tests. Then, the input of the functions to perform the 
#'   statistical tools is a matrix \code{XM_T} where each column corresponds to
#'   a vector formed by the values of a series \eqn{X_t}, from \eqn{t=1,\ldots,T},
#'   so that each row of the matrix correspond to a time \eqn{t}.
#'
#'   In  many real problems, such as those related to environmental phenomena, 
#'   the series of variables to analyze show a seasonal behaviour, and only one
#'   realization is available. In order to be able to apply the suggested tools
#'   to detect the existence of a  trend, the seasonal component has to be 
#'   removed and a sample of \eqn{M} uncorrelated series has to be obtained, 
#'   with a high enough value \eqn{M}. Those  problems can be solved, 
#'   preprocessing the data adequately. A wide set of tools to carry out a 
#'   preliminary analysis and to preprocess daily data with a yearly seasonal 
#'   pattern are implemented: \code{\link{series_double}}, 
#'   \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}} and \code{\link{series_untie}}.
#'
#'   There is also available a set of functions to characterize the occurrence 
#'   time, the value and other statistics related to the occurrence of records:
#'   \code{\link{I.record}}, \code{\link{L.record}}, \code{\link{M.record}},
#'   \code{\link{N.record}}, \code{\link{Nmean.record}} and \code{\link{records}}.
#'
#'   All the tests are based on the occurrence of records, but two families 
#'   are distinguished: the first is based on the number of records up to time
#'   t, \eqn{N_t} (\code{\link{foster.test}}, \code{\link{N_normal.test}}).
#'   The second family is based on the record indicator random variables,
#'   \eqn{I_t} and record times, \eqn{L_i} (\code{\link{L_lm.test}}, 
#'   \code{\link{L_lr.test}}, \code{\link{L_global.test}}, 
#'   \code{\link{P_chisq.test}}, \code{\link{P_exactPB.test}}, 
#'   \code{\link{P_regression.test}}). All of them are distribution-free tests 
#'   in time series for trend in location based on the null hypothesis that the
#'   record indicators are independent and the probabilies of record at time 
#'   \eqn{t} are \eqn{p_t=1/t}. 
#'
#'   Finally, there is another set of functions aiming to plot different 
#'   features related to records: \code{\link{foster.plot}},
#'   \code{\link{L.plot}}, \code{\link{N.plot}} and 
#'   \code{\link{P_regression.plot}}.
#'
#'   All the tests and plots can be applied to both upper and lower records, 
#'   positive or negative trends in location, and forward and backward series; 
#'   using the corresponding value in the argument \code{record}, \code{trend},
#'   or function \code{\link{series_rev}}, respectively.
#' @author Jorge Castillo-Mateo, Ana C. Cebrián
#' @docType package
#' @name RecordTest-package
NULL
