#' @title Regression test on record probabilities
#' @importFrom stats pf lm
#' @description This function performs a test based on a regression on the
#'   record probabilities \eqn{p_t} to study the hypothesis of the classical 
#'   record model.
#' @details The null hypothesis of this regression test is that in all the 
#'   vectors (columns in matrix \code{XM_T}), the probability of record at time
#'   \eqn{t} is \eqn{1/t}, so that \eqn{t p_t = 1}. Then, hypothesis 
#'   \eqn{H_0:\,p_t = 1/t}, \eqn{t=2, ..., T} is equivalent to 
#'   \eqn{H_0:\,\beta_0 = 1, \, \beta_1 = 0} where \eqn{\beta_0} and 
#'   \eqn{\beta_1} are the coefficients of the regression model 
#'   \eqn{t p_t=\beta_0 + \beta_1 t}. The  model has to be estimated by 
#'   weighted least squares since the response is heteroskedastic.
#'
#'   The F statistic is used to compare the regression model under the null 
#'   hypothesis and a linear regression model with no restriction (the 
#'   alterantive hypothesis is then that \eqn{t p_t} is a linear function of 
#'   time). This alternative hypothesis may be reasonable in many real 
#'   examples, but not always.
#'
#' @param XM_T A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of records to be 
#'   calculated, "upper" or "lower".
#' @return A \code{"htest"} object with elements:
#'   \item{null.value}{Value of \eqn{\beta_0} and \eqn{\beta_1} 
#'     under the null hypothesis.}
#'   \item{alternative}{Character string indicating the type of alternative
#'     hypothesis (two-sided).}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{estimate}{Value of \eqn{\hat \beta_0} and \eqn{\hat \beta_1}.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   \item{statistic}{Value of the F statistic.}
#'   \item{parameters}{Degrees of freedom of the F statistic.}
#'   \item{p.value}{P-value.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{P_exactPB.test}}, \code{\link{P_chisq.test}}, 
#'   \code{\link{P_regression.plot}}
#' @examples
#' P_regression.test(ZaragozaSeries)
#' @export P_regression.test
#'

P_regression.test <- function(XM_T, record = c('upper', 'lower')) {
  
  record <- match.arg(record)
  
  DNAME  <- deparse(substitute(XM_T))
  METHOD <- ifelse(record=='upper', "Regression test on upper record probabilities", 
                                    "Regression test on lower record probabilities")

  Trows <- NROW(XM_T)
  Mcols <- NCOL(XM_T)
  
  t <- 2:Trows

  # weigths
  VAR <- Mcols / (t - 1)
  P <- P.record(XM_T, record = record)[-1]
  t.P <- t * P

  model <- stats::lm(t.P ~ t, weights = VAR, y = TRUE)

  values <- linearHypothesis(model)
  
  structure(list(null.value = c("(intercept)" = 1, "t" = 0),
                 alternative = "two-sided",
                 method = METHOD,
                 estimate = model$coefficients,
                 data.name = DNAME,
                 statistic = values[2],
                 parameters = values[3:4],
                 p.value = values[1]), class = 'htest')
}

linearHypothesis <- function(lm) {
  
  df <- lm$df.residual
  
  RSS1 <- sum(lm$weights * lm$residuals^2)
  RSS0 <- sum(lm$weights * (1 - lm$y)^2)
  
  observed_value <- df * ( RSS0 - RSS1 ) / ( 2 * RSS1 )
  
  values <- c(stats::pf(observed_value, df1 = 2, df2 = df, lower.tail = FALSE), observed_value, 2, df)
  
  names(values) <- c('p-value', 'F', 'df1', 'df2')
  
  return(values)
}