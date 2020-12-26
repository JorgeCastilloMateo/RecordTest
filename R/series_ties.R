#' @title Summary of Record Ties
#' @description This function compares the number of strog and weak records
#'   to quantify whether rounding effects could greatly skew the conclusions.
#' @param X A numeric vector, matrix (or data frame).
#' @param record A character string indicating the type of record to be 
#'   assessed, "upper" or "lower".
#' @details This function is used in the data preparation (or pre-processing) 
#'   often required to apply the exploratory and inference tools based on 
#'   theory of records within this package.
#'  
#'  The theory of records on which the hypothesis tests are based assumes 
#'  that the random variables are continuous, proving that the probability 
#'  that two observations take the same value is zero. Most of the data 
#'  collected is rounded, giving a certain probability to the tie between 
#'  records, thereby reducing the number of new records(see, e.g., Wergen 
#'  et al. 2012). 
#'  
#'  This function summarises the difference between the number of observed 
#'  strong records and the weak records.
#'
#' @return A \code{list} object with elements:
#'   \item{`Number of records:`}{A vector containing the observed total, strong
#'     and weak number of records and the expected under IID.}
#'   \item{`\% of weak records:`}{Percentage of weak records within the total.}
#'   \item{`\% of weak records by position:`}{A vector with the percentage of
#'     weak records with \code{names} corresponding to its observed instant.}
#' @author Jorge Castillo-Mateo
#' @seealso \code{\link{series_double}}, \code{\link{series_record}}, 
#'   \code{\link{series_rev}}, \code{\link{series_split}}, 
#'   \code{\link{series_uncor}}, \code{\link{series_untie}} 
#' @references 
#' Wergen G, Volovik D, Redner S, Krug J (2012). 
#' “Rounding Effects in Record Statistics.”
#' \emph{Physical Review Letters}, \strong{109}(16), 164102. 
#' doi:\href{https://doi.org/10.1103/PhysRevLett.109.164102}{10.1103/PhysRevLett.109.164102}.
#' @examples
#' series_ties(ZaragozaSeries)
#' 
#' @export series_ties
series_ties <- function(X, record = c("upper", "lower")) {
  
  record <- match.arg(record)
  
  Trows <- NROW(X)
  Mcols <- NCOL(X)
  
  I_expected <- Mcols * sum(1 / 1:Trows)
  
  I_strong <- I.record(X, record = record, weak = FALSE)
  I_total  <- I.record(X, record = record, weak = TRUE)
  
  I_weak <- rowSums(I_total - I_strong)
  names(I_weak) <- 1:Trows

  sum_total <- c(sum(I_total), sum(I_strong), sum(I_weak),  I_expected)
  
  percent <- sum_total[3] / sum_total[1] * 100
  
  names(sum_total) <- c("Total", "Strong", "Weak", "Expected under IID")
  
  return(list("Number of records:" = sum_total, 
              "% of weak records:" = percent, 
              "% of weak records by position:" = I_weak / rowSums(I_total) * 100))
}