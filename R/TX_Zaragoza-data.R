#' @title Time Series of Daily Maximum Temperature at Zaragoza (Spain)
#' @description A dataset containing the series of daily maximum temperature at 
#'   Zaragoza aeropuerto (Spain), from 01/01/1951 to 31/12/2020. Zaragoza is 
#'   located at the north-east (+41:39:42 N, -001:00:29 W) of Iberian Peninsula
#'   at 247 m above mean sea level. This series is obtained from the ECA&D 
#'   series but it has been transformed, by removing days February
#'   29th and imputing the three missing values in the original series (by
#'   exponential weighted moving average expanded four days to each side). The
#'   variables are the following:
#'
#' \itemize{
#'   \item STAID : Station identifier
#'   \item SOUID : Source identifier
#'   \item DATE  : Date YYYYMMDD
#'   \item TX    : Maximum temperature in 0.1 ºC
#'   \item Q_TX  : quality code for TX (0='valid'; 1='suspect'; 9='missing')
#' }
#'
#' @docType data
#' @keywords datasets
#' @name TX_Zaragoza
#' @source \href{https://www.ecad.eu}{EUROPEAN CLIMATE ASSESSMENT & DATASET (ECA&D)}
#' @references 
#' Klein Tank AMG and Coauthors (2002). 
#' Daily Dataset of 20th-Century Surface Air Temperature and Precipitation 
#' Series for the European Climate Assessment.
#' \emph{International Journal of Climatology}, \strong{22}(12), 1441-1453. 
#' @usage data(TX_Zaragoza)
#' @seealso \code{\link{ZaragozaSeries}}
#' @format A data frame with 25550 rows and 5 variables.
"TX_Zaragoza"
