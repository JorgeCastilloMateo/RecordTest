#' @title Record Values and Record Times
#' 
#' @importFrom ggplot2 ggplot aes geom_step geom_point scale_colour_manual 
#'   scale_alpha_manual theme_bw theme labs
#'   
#' @description This function identifies (and plots if argument 
#'   \code{plot = TRUE}) the record values (\eqn{R_i}), and the record times 
#'   (\eqn{L_i}) in a vector, for all upper and lower records in forward and
#'   backward directions. 
#'   
#' @details 
#'   Customarily, the records in a time series \eqn{(X_t)} observed in \eqn{T} 
#'   instances \eqn{t = 1, 2, \ldots, T} can be obtained using chronological 
#'   order. Besides, we could also compute the records in similar sequences of 
#'   random variables if we consider reversed chronological order starting 
#'   from the last observation, i.e., \eqn{t' = T, \ldots, 2, 1}. The analysis 
#'   of series with reversed order is customarily referred to as backward, as 
#'   opposed to a forward analysis. 
#'   
#' @param X A numeric vector.
#' @param plot Logical. If \code{TRUE} (the default) the records are plotted.
#' @param direction A character string indicating the type of record to show 
#'   in the plot if \code{plot == TRUE}: \code{"forward"}, \code{"backward"} or
#'   \code{"both"} (see Details). 
#' @param variable Optional. A vector, containing other variable related 
#'   to \code{X} and measured at the same times. Only used if 
#'   \code{plot = FALSE}.
#' @param type Character string indicating if \code{X} is shown with 
#'   \code{"lines"} or \code{"points"}.
#' @param col,alpha Character and numeric vectors of length four, respectively.
#'   These arguments represent respectively the colour and transparency of the
#'   points or lines: trivial record, upper records, lower records and 
#'   observations respectively. 
#'   Vector names in the default are only indicative.
#' @param shape If \code{type == "points"}. 
#'   Integer vector of length 3 indicating the shape of the points
#'   for forward records, backward records and observations.
#'   Vector names in the default are only indicative.
#' @param linetype Integer vector of length 2 indicating the line type of the
#'   step functions in the forward and backward records, respectively.
#'   Vector names in the default are only indicative.
#' @return If \code{plot = TRUE} a ggplot object, otherwise a list with four 
#'   data frames where the first column are the record times, the second the 
#'   record values and, if \code{variable} is not null, the third column are 
#'   their values at the record times, respectively for upper and lower records
#'   in forward and backward series.
#'   
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{I.record}}, \code{\link{series_double}}, 
#'   \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}}
#' @examples
#' Y <- c(5, 7, 3, 6, 19, 2, 20)
#' records(Y, plot = FALSE, variable = seq_along(Y))
#' 
#' # Show the whole series and its upper and lower records
#' records(TX_Zaragoza$TX)
#' # Compute tables for the whole series
#' TxZ.record <- records(TX_Zaragoza$TX, plot = FALSE, variable = TX_Zaragoza$DATE)
#' TxZ.record
#' names(TxZ.record)
#' # To show the Forward Upper records
#' TxZ.record[[1]]
#' plot(TxZ.record[[1]]$Times, TxZ.record[[1]]$Values)
#' 
#' # Annual maximum daily maximum temperatures
#' TxZ <- apply(series_split(TX_Zaragoza$TX), 1, max)
#' # Plot for the records in forward and backward directions
#' records(TxZ, direction = "both")
#' # Compute tables for the annual maximum 
#' records(TxZ, plot = FALSE, variable = 1951:2020)
#' 
#' @export records
records <- function(X, 
                    plot = TRUE, 
                    direction = c("forward", "backward", "both"),
                    variable, 
                    type = c("lines", "points"),
                    col = c("T" = "black", "U" = "salmon", "L" = "skyblue", "O" = "black"), 
                    alpha = c("T" = 1, "U" = 1, "L" = 1, "O" = 1),
                    shape = c("F" = 19, "B" = 4, "O" = 19),
                    linetype = c("F" = 1, "B" = 2)) {
  
  direction <- match.arg(direction)
  type <- match.arg(type)
  
  names(col)   <- NULL
  names(alpha) <- NULL
  names(shape) <- NULL
  
  Tlength <- length(X)
  
  Xrev <- rev(X)
  I.FU <- I.record(X, record = "upper")
  I.FL <- I.record(X, record = "lower")
  I.BU <- I.record(Xrev, record = "upper")
  I.BL <- I.record(Xrev, record = "lower")
  
  if (plot) {
    
    # ggplot
    graf <- ggplot2::ggplot(data.frame(X), mapping = ggplot2::aes(x = seq_len(Tlength))) + 
      ggplot2::scale_colour_manual(name = "Records", values = c("Forward Trivial" = col[1], "Forward Upper" = col[2], "Observation" = col[4], "Forward Lower" = col[3],
                                                               "Backward Trivial" = col[1], "Backward Upper" = col[2], "Observation" = col[4], "Backward Lower" = col[3], "NA" = 1), 
                                   breaks = c("Forward Trivial", "Backward Trivial", "Forward Upper", "Backward Upper", "Forward Lower", "Backward Lower", "Observation")) +
      ggplot2::scale_alpha_manual(name = "Records", values = c("Forward Trivial" = alpha[1], "Forward Upper" = alpha[2], "Observation" = alpha[4], "Forward Lower" = alpha[3],
                                                              "Backward Trivial" = alpha[1], "Backward Upper" = alpha[2], "Observation" = alpha[4], "Backward Lower" = alpha[3], "NA" = 0), 
                                  breaks = c("Forward Trivial", "Backward Trivial", "Forward Upper", "Backward Upper", "Forward Lower", "Backward Lower", "Observation")) +
      ggplot2::scale_shape_manual(name = "Records", values = c("Forward Trivial" = shape[1], "Forward Upper" = shape[1], "Observation" = shape[3], "Forward Lower" = shape[1],
                                                              "Backward Trivial" = shape[2], "Backward Upper" = shape[2], "Observation" = shape[3], "Backward Lower" = shape[2], "NA" = 1), 
                                  breaks = c("Forward Trivial", "Backward Trivial", "Forward Upper", "Backward Upper", "Forward Lower", "Backward Lower", "Observation")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(title = "Record values and record times", 
                    subtitle = paste("Data:", deparse(substitute(X))), 
                    x = "t", y = expression(X[t]))
    
    if (direction == "forward") {
      # Forward
      I.F <- I.FU - I.FL
      recordFactor <- as.factor(I.F)
      levels(recordFactor) <- c(levels(recordFactor), "Forward Trivial")
      recordFactor[1] <- "Forward Trivial"
      levels(recordFactor)[levels(recordFactor) == "-1"] <- "Forward Lower"
      levels(recordFactor)[levels(recordFactor) == "0"]  <- "Observation"
      levels(recordFactor)[levels(recordFactor) == "1"]  <- "Forward Upper"
      Xna <- X
      Xna[is.na(X)] <- X[1]
      recordUpper <- cummax(Xna)
      recordLower <- cummin(Xna)
      
      graf <- graf +
        ggplot2::geom_step(ggplot2::aes(y = recordUpper), direction = "hv", colour = col[2], alpha = alpha[2], linetype = linetype[1]) +
        ggplot2::geom_step(ggplot2::aes(y = recordLower), direction = "hv", colour = col[3], alpha = alpha[3], linetype = linetype[1])
      if (type == "lines") {
        graf <- graf + 
          ggplot2::geom_line(ggplot2::aes(y = X), colour = col[4], alpha = alpha[4])
      } else { # type == "points
        graf <- graf + 
          ggplot2::geom_point(ggplot2::aes(y = X, colour = recordFactor, alpha = recordFactor, shape = recordFactor))
      }
        
    } else if (direction == "backward") {
      # Backward
      I.B <- I.BU - I.BL
      recordFactor.B <- as.factor(I.B)
      levels(recordFactor.B) <- c(levels(recordFactor.B), "Backward Trivial")
      recordFactor.B[1] <- "Backward Trivial"
      levels(recordFactor.B)[levels(recordFactor.B) == "-1"] <- "Backward Lower"
      levels(recordFactor.B)[levels(recordFactor.B) == "0"]  <- "Observation"
      levels(recordFactor.B)[levels(recordFactor.B) == "1"]  <- "Backward Upper"
      Xrevna <- Xrev
      Xrevna[is.na(Xrevna)] <- Xrev[1]
      recordUpper.B <- cummax(Xrevna)
      recordLower.B <- cummin(Xrevna)
      
      graf <- graf +
        ggplot2::geom_step(ggplot2::aes(x = rev(seq_len(Tlength)), y = recordUpper.B), direction = "vh", colour = col[2], alpha = alpha[2], linetype = linetype[2]) +
        ggplot2::geom_step(ggplot2::aes(x = rev(seq_len(Tlength)), y = recordLower.B), direction = "vh", colour = col[3], alpha = alpha[3], linetype = linetype[2])
      if (type == "lines") {
        graf <- graf + 
          ggplot2::geom_line(ggplot2::aes(x = rev(seq_len(Tlength)), y = rev(X)), colour = col[4], alpha = alpha[4])
      } else { # type == "points
        graf <- graf + 
          ggplot2::geom_point(ggplot2::aes(x = rev(seq_len(Tlength)), y = rev(X), colour = recordFactor.B, alpha = recordFactor.B, shape = recordFactor.B)) 
      }
      
    } else { # direction == "both"
      # Forward
      I.F <- I.FU - I.FL
      recordFactor <- as.factor(I.F)
      levels(recordFactor) <- c(levels(recordFactor), "Forward Trivial", "NA")
      recordFactor[1] <- "Forward Trivial"
      levels(recordFactor)[levels(recordFactor) == "-1"] <- "Forward Lower"
      levels(recordFactor)[levels(recordFactor) == "0"]  <- "Observation"
      levels(recordFactor)[levels(recordFactor) == "1"]  <- "Forward Upper"
      Xna <- X
      Xna[is.na(Xna)] <- X[1]
      recordUpper <- cummax(Xna)
      recordLower <- cummin(Xna)
      # Backward
      I.B <- I.BU - I.BL
      recordFactor.B <- as.factor(I.B)
      levels(recordFactor.B) <- c(levels(recordFactor.B), "Backward Trivial", "NA")
      recordFactor.B[1] <- "Backward Trivial"
      levels(recordFactor.B)[levels(recordFactor.B) == "-1"] <- "Backward Lower"
      levels(recordFactor.B)[levels(recordFactor.B) == "0"]  <- "Observation"
      levels(recordFactor.B)[levels(recordFactor.B) == "1"]  <- "Backward Upper"
      Xrevna <- Xrev
      Xrevna[is.na(Xrevna)] <- Xrev[1]
      recordUpper.B <- cummax(Xrevna)
      recordLower.B <- cummin(Xrevna)
      remove.F <- as.logical((recordFactor == "Observation") + (rev(recordFactor.B) != "Observation") == 2)
      remove.B <- as.logical((rev(recordFactor) != "Observation") + (recordFactor.B == "Observation") == 2)
      recordFactor[remove.F] <- "NA"
      recordFactor.B[remove.B] <- "NA"

      graf <- graf +
        ggplot2::geom_step(ggplot2::aes(x = rev(seq_len(Tlength)), y = recordUpper.B), direction = "vh", colour = col[2], alpha = alpha[2], linetype = linetype[2]) +
        ggplot2::geom_step(ggplot2::aes(x = rev(seq_len(Tlength)), y = recordLower.B), direction = "vh", colour = col[3], alpha = alpha[3], linetype = linetype[2]) +
        ggplot2::geom_step(ggplot2::aes(y = recordUpper), direction = "hv", colour = col[2], alpha = alpha[2], linetype = linetype[1]) +
        ggplot2::geom_step(ggplot2::aes(y = recordLower), direction = "hv", colour = col[3], alpha = alpha[3], linetype = linetype[1])
      if (type == "lines") {
        graf <- graf + 
          ggplot2::geom_line(ggplot2::aes(y = X), colour = col[4], alpha = alpha[4])
      } else { # type == "points
        graf <- graf + 
          ggplot2::geom_point(ggplot2::aes(x = rev(seq_len(Tlength)), y = rev(X), colour = recordFactor.B, alpha = recordFactor.B, shape = recordFactor.B)) +
          ggplot2::geom_point(ggplot2::aes(y = X, colour = recordFactor, alpha = recordFactor, shape = recordFactor))
      }
    }
    
    return(graf)
    
  } else {
    
    if (!missing(variable) && (Tlength != length(variable))) {
      
      stop("'X' and 'variable' should have the same length")
    }
    
    L.FU <- which(I.FU == 1)
    R.FU <- X[L.FU]
    L.FL <- which(I.FL == 1)
    R.FL <- X[L.FL]
    L.BU <- which(I.BU == 1)
    R.BU <- rev(X)[L.BU]
    L.BL <- which(I.BL == 1)
    R.BL <- rev(X)[L.BL]
  }
  
  if (!missing(variable)) {
    
    variable <- as.character(variable)
    Var_L.FU <- variable[L.FU]
    Var_L.FL <- variable[L.FL]
    Var_L.BU <- variable[Tlength - L.BU + 1]
    Var_L.BL <- variable[Tlength - L.BL + 1]
    table <- list("Forward.Upper.record"  = data.frame("Times" = L.FU, "Values" = R.FU, "Variable" = Var_L.FU), 
                  "Forward.Lower.record"  = data.frame("Times" = L.FL, "Values" = R.FL, "Variable" = Var_L.FL),
                  "Backward.Upper.record" = data.frame("Times" = L.BU, "Values" = R.BU, "Variable" = Var_L.BU), 
                  "Backward.Lower.record" = data.frame("Times" = L.BL, "Values" = R.BL, "Variable" = Var_L.BL))
    
  } else {
    
    table <- list("Forward.Upper.record"  = data.frame("Times" = L.FU, "Values" = R.FU), 
                  "Forward.Lower.record"  = data.frame("Times" = L.FL, "Values" = R.FL),
                  "Backward.Upper.record" = data.frame("Times" = L.BU, "Values" = R.BU), 
                  "Backward.Lower.record" = data.frame("Times" = L.BL, "Values" = R.BL))
  }
  
  return(table)
}
