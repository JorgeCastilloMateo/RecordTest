#' @title Regression Test for Record Probabilities
#' @importFrom stats pf lm
#' @description This function performs a linear hypothesis test based on a
#'   regression for the record probabilities \eqn{p_t} to study the hypothesis
#'   of the classical record model.
#' @details 
#'   The null hypothesis is that the data come from a population with 
#'   independent and identically distributed realisations. This implies that
#'   in all the vectors (columns in matrix \code{X}), the sample probability 
#'   of record at time \eqn{t} (\code{\link{p.record}}) is \eqn{1/t}, so that
#'   \deqn{t \, \textrm{E}(\hat p_t) = 1.} 
#'   Then, 
#'   \deqn{H_0:\,p_t = 1/t, \, t=2, ..., T \iff H_0:\,\beta_0 = 1, \, \beta_1 = 0,} 
#'   where \eqn{\beta_0} and \eqn{\beta_1} are the coefficients of the 
#'   regression model 
#'   \deqn{t \, \textrm{E}(\hat p_t) = \beta_0 + \beta_1 t.} 
#'   The  model has to be estimated by weighted least squares since the 
#'   response is heteroskedastic.
#'   
#'   Other models can be considered with the \code{formula} argument. 
#'   However, for the test to be correct, the model that assigns 1 to all 
#'   responses must be nested in the bigger one, either leaving the intercept 
#'   free or setting the intercept to 1 (see Examples for possible models).
#'
#'   The \eqn{F} statistic is computed for carrying out a Wald-test-based 
#'   comparison between the restricted model under the null hypothesis and
#'   the more general model (e.g., the alterantive hypothesis 
#'   where \eqn{t \, \textrm{E}(\hat p_t)} is a linear function of time \eqn{t}). 
#'   This alternative hypothesis may be reasonable in many real examples, 
#'   but not always.
#'   
#'   If the sample size (i.e., the number of series or columns of \code{X})
#'   is lower than 8 or 12 the distribution \eqn{F} is not fulfilled, so the 
#'   \code{simulate.p.value} option is recommended in this case.
#'
#' @param X A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of records to be 
#'   calculated, "upper" or "lower".
#' @param formula "\code{\link{formula}}" to use in \code{\link{lm}} function, 
#'   e.g., \code{y ~ x}, \code{y ~ poly(x, 2, raw = TRUE)}, \code{y ~ log(x)}.
#'   By default \code{formula = y ~ x}. See Note for a caveat.
#' @param simulate.p.value Logical. Indicates whether to compute p-values by
#'   Monte Carlo simulation. It is recommended if the number of columns of
#'   \code{X} (i.e., the number of series) is lower than 12, since for lower 
#'   values the size of the test is not fulfilled.    
#' @param B If \code{simulate.p.value = TRUE}, an integer specifying the 
#'   number of replicates used in the Monte Carlo estimation.
#' @return A \code{"htest"} object with elements:
#'   \item{null.value}{Value of the coefficients under the null hypothesis
#'     when more than one coefficient is fitted.}
#'   \item{alternative}{Character string indicating the type of alternative
#'     hypothesis.}
#'   \item{method}{A character string indicating the type of test performed.}
#'   \item{estimate}{Value of the fitted coefficients.}
#'   \item{data.name}{A character string giving the name of the data.}
#'   \item{statistic}{Value of the \eqn{F} statistic.}
#'   \item{parameters}{Degrees of freedom of the \eqn{F} statistic.}
#'   \item{p.value}{P-value.}
#' @note IMPORTANT: In \code{formula} the intercept has to be free or fixed
#'   to 1 so that the test is correct.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{chisq.test}}, \code{\link{p.plot}}
#' @examples
#' # Simple test for upper records (p-value = 0.01047)
#' p.test(ZaragozaSeries)
#' # Simple test for lower records (p-value = 9.178e-05)
#' p.test(ZaragozaSeries, record = "lower")
#' 
#' # Fit a 2nd term polynomial for upper records (p-value = 0.01187)
#' p.test(ZaragozaSeries, formula = y ~ I(x^2))
#' # Fit a 2nd term polynomial for lower records (p-value = 8.007e-05)
#' p.test(ZaragozaSeries, record = "lower", formula = y ~ I(x^2))
#' 
#' # Fix the intercept to 1 for upper records (p-value = 0.005557)
#' p.test(ZaragozaSeries, formula = y ~ I(x-1) - 1 + offset(rep(1, length(x))))
#' # Fix the intercept to 1 for lower records (p-value = 2.467e-05)
#' p.test(ZaragozaSeries, record = "lower", formula = y ~ I(x-1) - 1 + offset(rep(1, length(x))))
#' 
#' # Simulate p-value when the number of series is small
#' TxZ <- apply(series_split(TX_Zaragoza$TX), 1, max)
#' p.test(TxZ, simulate.p.value = TRUE)
#' @export p.test

p.test <- function(X, 
                   record = c("upper", "lower"), 
                   formula = y ~ x, 
                   simulate.p.value = FALSE, 
                   B = 1000) {
  
  record <- match.arg(record)
  
  DNAME  <- deparse(substitute(X))
  METHOD <- paste("Regression test on", record, "record probabilities")
  
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  
  t <- 2:Trows
  
  # weigths
  invVAR <- Mcols / (t - 1)
  p <- p.record(X, record = record)[-1]
  tp <- t * p

  model <- stats::lm(data = data.frame(x = t, y = tp, invVAR), formula = formula, weights = invVAR, y = TRUE)
  values <- .linearHypothesis(model)
  
  p.fun <- function(y) {
    model <- stats::lm(data = data.frame(x = t, y, invVAR), formula = formula, weights = invVAR, y = TRUE)
    values <- .linearHypothesis(model)[1]
    return(values)
  }
  
  if (simulate.p.value) {
    METHOD <- paste(METHOD, "with simulated p-value (based on", B, "replicates)")
    tpB <- t * matrix(stats::rbinom(n = (Trows - 1) * B, size = Mcols, prob = 1 / t), ncol = B) / Mcols
    FstatB <- apply(tpB, 2, p.fun)
    pv <- sum(FstatB >= values[1]) / B
  } else {
    pv <- stats::pf(values[1], df1 = values[2], df2 = values[3], lower.tail = FALSE)
  }
  
  H0 <- model$coefficients * 0
  
  if (length(H0) == 1) {
    alternative <- ifelse(names(H0) == "(Intercept)",
                          "true (Intercept) is not 1",
                          paste("true", names(H0), "is not 0"))
    structure(list(alternative = alternative,
                   method = METHOD,
                   estimate = model$coefficients,
                   data.name = DNAME,
                   statistic = values[1],
                   parameters = values[2:3],
                   p.value = pv), class = "htest")
  } else {
    H0["(Intercept)"] <- 1
    structure(list(null.value = H0,
                   alternative = "two-sided for record probabilities",
                   method = METHOD,
                   estimate = model$coefficients,
                   data.name = DNAME,
                   statistic = values[1],
                   parameters = values[2:3],
                   p.value = pv), class = "htest")
  }
}

.linearHypothesis <- function(lm) {
  
  df <- lm$df.residual
  diff_df <- lm$rank
  
  RSS1 <- sum(lm$weights * lm$residuals^2)
  RSS0 <- sum(lm$weights * (1 - lm$y)^2)
  
  observed_value <- df * ( RSS0 - RSS1 ) / ( diff_df * RSS1 )
  
  values <- c(observed_value, diff_df, df)
  
  names(values) <- c("F", "df1", "df2")
  
  return(values)
}