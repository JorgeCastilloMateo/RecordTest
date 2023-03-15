#' @title \strong{RecordTest}: A Package for Testing the Classical Record Model
#' @aliases RecordTest-package
#' @aliases RecordTest
#' @description \strong{RecordTest} provides data preparation, exploratory data 
#'   analysis and inference tools based on theory of records to describe the 
#'   record occurrence and detect trends, change-points or non-stationarities 
#'   in the tails of the time series. Details about the implemented tools can
#'   be found in Castillo-Mateo, Cebrián and Asín (2023).
#' @details 
#'   The Classical Record Model:
#' 
#'   Record statistics are used primarily to quantify the stochastic behaviour
#'   of a process at never-seen-before values, either upper or lower. The setup
#'   of independent and identically distibuted (IID) continuous random 
#'   variables (RVs), often called the classical record model, is 
#'   particularly interesting because the common continuous distribution 
#'   underlying the IID continuous RVs will not affect the distribution of the variables 
#'   relative to the record occurrence.  
#'   Many fields have begun to use the theory of records to study these 
#'   remarkable events. Particularly productive is the study of 
#'   record-breaking temperatures and their connection with climate change, 
#'   but also records in other environmental fields (precipitations, floods, 
#'   earthquakes, etc.), economy, biology, physics or even sports have been 
#'   analysed.
#'   See Arnold, Balakrishnan and Nagaraja (1998) for an extensive theoretical
#'   introduction to the theory of records and in particular the classical 
#'   record model. See Foster and Stuart (1954), Diersen and Trenkler (1996, 
#'   2001) and Cebrián, Castillo-Mateo and Asín (2022) for the 
#'   distribution-free trend detection tests and Castillo-Mateo (2022) for the
#'   distribution-free change-point detection tests based on the classical 
#'   record model. For an easy introduction to \strong{RecordTest} use 
#'   \code{vignette("RecordTest")}.
#' 
#'   This package provides tests to study the hypothesis of the classical 
#'   record model, that is that the record occurrence from a series of values 
#'   observed at regular time units come from an IID series of continuous RVs. 
#'   If we have sequences of independent variables with no seasonal component,
#'   the hypothesis of IID variables is equivalent to test the hypothesis of 
#'   homogeneity and stationarity.
#' 
#'   The functions in the data preparation step:
#'   
#'   The functions admit a vector \code{X} corresponding to a single series as
#'   an argument. However, some situations could take advantage of having 
#'   \eqn{M} uncorrelated vectors to infer from the sample. Then, the input of
#'   the functions to perform the statistical tools can be a matrix \code{X} 
#'   where each column corresponds to a vector formed by the values of a 
#'   series \eqn{X_t}, for \eqn{t=1,\ldots,T}, so that each row of the matrix
#'   correspond to a time \eqn{t}.
#'
#'   In  many real problems, such as those related to environmental phenomena, 
#'   the series of variables to analyse show a seasonal behaviour, and only one
#'   realisation is available. In order to be able to apply the suggested tools
#'   to detect the existence of a trend, the seasonal component has to be 
#'   removed and a sample of \eqn{M} uncorrelated series should be obtained. 
#'   Those problems can be solved by preparing the data adequately. 
#'   A wide set of tools to carry out a preliminary analysis and to prepare 
#'   data with a seasonal pattern are implemented in the following functions.
#'   
#'   \code{\link{series_record}}: If only the record times are available.
#'   
#'   \code{\link{series_split}}, \code{\link{series_double}}: To split the
#'   series in several subseries and remove the seasonal component and 
#'   autocorrelation.
#'   
#'   \code{\link{series_uncor}}: To extract a subset of uncorrelated subseries 
#'   
#'   \code{\link{series_ties}}, \code{\link{series_untie}}: To deal with record
#'   ties.
#'   
#'   \code{\link{series_rev}}: To study the series backwards.
#'
#'   The functions to compute the record statistics are:
#'   
#'   \code{\link{I.record}}: Computes the observed record indicators. \code{NA} 
#'   values are taken as no records unless they appear at \eqn{t = 1}.
#'   
#'   \code{\link{N.record}}, \code{\link{Nmean.record}}: Compute the observed
#'   number of records up to time \eqn{t}.
#'   
#'   \code{\link{S.record}}: Computes the observed number of records at every
#'   time \eqn{t}, using \eqn{M} series.
#'   
#'   \code{\link{p.record}}: Computes the estimated record probability at every
#'   time \eqn{t}, using \eqn{M} series.
#'   
#'   \code{\link{L.record}}: Computes the observed record times.
#'   
#'   \code{\link{R.record}}: Computes the observed record values.
#'   
#'   The functions to compute the tests:
#'   
#'   All the tests performed are distribution-free/non-parametric tests in 
#'   time series for trend, change-point and non-stationarity in the extremes 
#'   of the distribution based on the null hypothesis that the record 
#'   indicators are independent and the probabilities of record at time \eqn{t}
#'   are \eqn{p_t = 1 / t}. 
#'   
#'   \code{\link{change.point}}: Implements Castillo-Mateo change-point tests.
#'   
#'   \code{\link{foster.test}}: Implements Foster-Stuart and Diersen-Trenkler
#'   trend tests. 
#'   
#'   \code{\link{N.test}}: Implements tests based on the (weighted) number of
#'   records.
#'   
#'   \code{\link{brown.method}}: Brown's method to combine dependent p-values 
#'   from \code{\link{N.test}}.
#'   
#'   \code{\link{fisher.method}}: General function to apply Fisher's method to
#'   independent p-values.
#'   
#'   \code{\link{p.regression.test}}: Implements a regression test based on the
#'   record probabilities.
#'   
#'   \code{\link{p.chisq.test}}: Implements a \eqn{\chi^2}-test based on the 
#'   record probabilities.
#'   
#'   \code{\link{lr.test}}: Implements likelihood ratio tests based on the 
#'   record indicators.
#'   
#'   \code{\link{score.test}}: Implements score or Lagrange multiplier
#'   tests based on the record indicators.
#'
#'   The functions to compute the graphical tools:
#'   
#'   \code{\link{records}}: Shows the series remarking its records.
#'   
#'   \code{\link{L.plot}}: Shows record times in several series.
#'   
#'   \code{\link{foster.plot}}: Shows plots based on Foster-Stuart and
#'   Diersen-Trenkler statistics.
#'   
#'   \code{\link{N.plot}}: Shows the (weighted) number of records.
#'    
#'   \code{\link{p.plot}}: Shows the record probabilities in different plots.
#'
#'   All the tests and graphical tools can be applied to both upper and lower 
#'   records in the forward and backward directions.
#'   
#'   Other functions:
#'   
#'   \code{\link{rcrm}}: Random generation for the classical record model.
#'   
#'   \code{\link{dpoisbinom}}, \code{\link{ppoisbinom}}, 
#'   \code{\link{qpoisbinom}}, \code{\link{rpoisbinom}}: Density, distribution
#'   function, quantile function and random generation for the Poisson binomial
#'   distribution. Related to the probability distribution function of the 
#'   number of records under the null hypothesis.
#'   
#'   Example datasets:
#'
#'   There are two example datasets included with this package. It is possible
#'   to load these datasets into R using the \code{data} function. The 
#'   datasets have their own help file, which can be accessed by 
#'   \code{help([dataset_name])}. 
#'   Data included with \strong{RecordTest} are:
#'   
#'   \code{\link{TX_Zaragoza}} - Daily maximum temperatures at Zaragoza 
#'   (Spain).
#'   
#'   \code{\link{ZaragozaSeries}} - Split and uncorrelated subseries 
#'   \code{\link{TX_Zaragoza}$TX}.
#'   
#'   \code{\link{Olympic_records_200m}} - 200-meter Olympic records from 1900
#'   to 2020.
#'   
#'   To see how to cite \strong{RecordTest} in publications or elsewhere,
#'   use \code{citation("RecordTest")}. 
#' @author Jorge Castillo-Mateo  <jorgecm@unizar.es>, AC Cebrián, J Asín
#' @references 
#' Arnold BC, Balakrishnan N, Nagaraja HN (1998). 
#' \emph{Records}. 
#' Wiley Series in Probability and Statistics. Wiley, New York.
#' 
#' Castillo-Mateo J (2022).
#' “Distribution-Free Changepoint Detection Tests Based on the Breaking of Records.”
#' \emph{Environmental and Ecological Statistics}, \strong{29}(3), 655-676. 
#' \doi{10.1007/s10651-022-00539-2}
#' 
#' Castillo-Mateo J, Cebrián AC, Asín J (2023).
#' “\strong{RecordTest}: An \code{R} Package to Analyze Non-Stationarity in the Extremes Based on Record-Breaking Events.”
#' \emph{Journal of Statistical Software}, \strong{106}(5), 1-28. 
#' \doi{10.18637/jss.v106.i05}
#'
#' Cebrián AC, Castillo-Mateo J, Asín J (2022).
#' “Record Tests to Detect Non Stationarity in the Tails with an Application to Climate Change.”
#' \emph{Stochastic Environmental Research and Risk Assessment}, \strong{36}(2): 313-330. 
#' \doi{10.1007/s00477-021-02122-w}
#' 
#' Diersen J, Trenkler G (1996). “Records Tests for Trend in Location.”
#' \emph{Statistics}, \strong{28}(1), 1-12.
#' \doi{10.1080/02331889708802543}
#' 
#' Diersen J, Trenkler G (2001). 
#' “Weighted Records Tests for Splitted Series of Observations.”
#' In J Kunert, G Trenkler (eds.), 
#' \emph{Mathematical Statistics with Applications in Biometry: Festschrift in Honour of Prof. Dr. Siegfried Schach}, 
#' pp. 163–178. Lohmar: Josef Eul Verlag.
#' 
#' Foster FG, Stuart A (1954). 
#' “Distribution-Free Tests in Time-Series Based on the Breaking of Records.”
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 
#' \strong{16}(1), 1-22.
#' @docType package
#' @name RecordTest-package
NULL
