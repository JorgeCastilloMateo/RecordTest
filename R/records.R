#' @title Record values and record times
#' @importFrom ggplot2 ggplot aes geom_step geom_point scale_colour_manual 
#'   scale_alpha_manual theme_bw theme labs
#' @description This function identifies (and plot if argument 
#'   \code{plot = TRUE}) the record values (\eqn{R_i}), and the record times 
#'   (\eqn{L_i}) in a vector, for both upper and lower records. 
#' @param X_T A numeric vector.
#' @param plot logical. If \code{TRUE} (the default) the records are plotted.
#' @param variables Optional. A matrix, containing other variables related 
#'   to \code{X_T} and measured at the same times. Only used if 
#'   \code{plot = FALSE}.
#' @param colour,alpha Character and numeric vectors of length 4, respectively.
#'   This arguments represent respectively the colour and transparency of the
#'   points: trivial record, upper records, observations and lower records, 
#'   respectively.
#' @return If \code{plot = TRUE} a ggplot object, otherwise a list with two 
#'   data frames where the first column are the record times, the second the 
#'   record values and, if \code{variables} is not null, the third column are 
#'   their values at the record times, respectively for upper and lower records.
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{series_double}}, 
#'   \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' Y <- c(5, 7, 3, 6, 19, 2, 20)
#' records(Y, plot = FALSE)
#' 
#' # compute tables
#' records(TX_Zaragoza$TX, plot = FALSE, variables = TX_Zaragoza$DATE)
#' # compute plot
#' records(TX_Zaragoza$TX, alpha = c(1, 1, 0.01, 1))
#' 
#' @export records

records <- function(X_T, plot = TRUE, variables = NULL, 
                    colour = c('black', 'skyblue2', 'limegreen', 'salmon'), 
                    alpha = c(1, 1, 0.1, 1)) {
  
  Tlength <- length(X_T)
  
  IU <- I.record(X_T, record = 'upper')
  IL <- I.record(X_T, record = 'lower')
  
  if (plot) {
    
    I <- IU - IL
    
    recordFactor <- as.factor(I)
    levels(recordFactor) <- c(levels(recordFactor), 'Trivial record')
    recordFactor[1] <- 'Trivial record'
    levels(recordFactor)[levels(recordFactor) == '-1'] <- 'Lower record'
    levels(recordFactor)[levels(recordFactor) == '0']  <- 'Observation'
    levels(recordFactor)[levels(recordFactor) == '1']  <- 'Upper record'
    
    recordUpper <- cummax(X_T)
    recordLower <- cummin(X_T)
    
    graf <- ggplot2::ggplot(data.frame(X_T), mapping = ggplot2::aes(x = seq_len(Tlength))) + 
      ggplot2::geom_step(ggplot2::aes(y = recordUpper), direction = "hv", colour = colour[2], alpha = 0.6) +
      ggplot2::geom_step(ggplot2::aes(y = recordLower), direction = "hv", colour = colour[4], alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(y = X_T, colour = recordFactor, alpha = recordFactor)) + 
      ggplot2::scale_colour_manual(name = "", values = c('Trivial record' = colour[1], 'Upper record' = colour[2], 'Observation' = colour[3], 'Lower record' = colour[4])) +
      ggplot2::scale_alpha_manual(name = "", values = c('Trivial record' = alpha[1], 'Upper record' = alpha[2], 'Observation' = alpha[3], 'Lower record' = alpha[4])) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'bottom') +
      ggplot2::labs(title = 'Record values and record times', 
                    subtitle = paste('Data:', deparse(substitute(X_T))), 
                    x = "t", y = expression(X[t]))
    
    return(graf)
    
  } else {
    
    if (!is.null(variables) && (Tlength != length(variables))) {
      
      stop("X_T and variables must have the same number of rows")
    }
    
    LU <- which(IU == 1)
    RU <- X_T[LU]
    
    LL <- which(IL == 1)
    RL <- X_T[LL]
  }
  
  if (!is.null(variables)) {
    
    variables <- as.character(variables)
    Var_LU <- variables[LU]
    Var_LL <- variables[LL]
    table <- list('Upper record' = data.frame('Times' = LU, 'Values' = RU, Var_LU), 
                  'Lower record' = data.frame('Times' = LL, 'Values' = RL, Var_LL))
    
  } else {
    
    table <- list('Upper record' = data.frame('Times' = LU, 'Values' = RU), 
                  'Lower record' = data.frame('Times' = LL, 'Values' = RL))
  }
  
  return(table)
}
