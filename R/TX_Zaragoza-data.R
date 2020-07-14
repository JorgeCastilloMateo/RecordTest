#' @title Time series of daily maximum temperature at Zaragoza (Spain)
#' @description A dataset containing the series of daily maximum temperature at 
#'   Zaragoza (Spain), from 01/01/1953 to 31/12/2018. This series is obtaned 
#'   from the ECA series but it has been transformed, by removing days February 
#'   29th and filling the missing values. The variables are the following:
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
#' @source \href{http://www.ecad.eu}{EUROPEAN CLIMATE ASSESSMENT & DATASET (ECA&D)}
#' @references Klein Tank, A.M.G. and Coauthors, 2002. 
#'   Daily dataset of 20th-century surface air temperature and precipitation 
#'   series for the European Climate Assessment.
#'   \emph{Int. J. of Climatol.}, 22, 1441-1453.
#' @usage data(TX_Zaragoza)
#' 
#' @format A data frame with 28670 rows and 5 variables.
"TX_Zaragoza"
