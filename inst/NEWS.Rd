\name{NEWS}
\title{News for Package \pkg{RecordTest}}
\newcommand{\cpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\issue}{\href{https://github.com/JorgeCastilloMateo/RecordTest/issues/#1}{(issue #1)}}

\section{Changes in Version 2.2.0 (?? ?? 2023)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Implementation of "approximate" or "random" permutation tests with the argument \code{permutation.test} in the functions \code{change.point}, \code{foster.test}, and \code{N.test}. The tests use samples with repetition of permutations, they do not use all possible permutations. Unlike previous exact, asymptotic, or Monte Carlo tests, permutation tests still ensure test size by introducing a matrix \code{X} with dependent columns as an argument.
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item \pkg{RecordTest} logo update.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Function \code{records} when \code{plot = FALSE} returned the wrong backward record values. It returned the values of the forward series on the backward record times. This issue was fixed.
      \item Function \code{L.record} returned a matrix instead of a list when all series had the same number of records. Function \code{R.record} returned incorrect output when all series had the same number of records or there were only one series. These issues were fixed and both functions return a list in any case.
      \item Function \code{series_untie} now admits \code{NA} or \code{NaN} values without returning \code{NaN}'s.
    }
  }
}

\section{Changes in Version 2.1.2 (15 March 2023)}{
  \subsection{MODIFICATIONS}{
    \itemize{
      \item \sQuote{CITATION} and references update for the new JSS publication of Castillo-Mateo, Cebrian and Asin (2023) <doi:10.18637/jss.v106.i05>.
    }
  }
}

\section{Changes in Version 2.1.1 (03 November 2022)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New dataset \code{TX_Zaragoza29F} that includes Zaragoza series with the observations from the 29th of February.
      \item Function \code{series_uncor} now admits any test introduced as an argument to test correlation or dependence across columns. Default is \code{cor.test}, which is the test used in previous versions.
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item References update for the publication of Cebrian, Castillo-Mateo and Asin (2022) <doi:10.1007/s00477-021-02122-w> and Castillo-Mateo (2022) <doi:10.1007/s10651-022-00539-2>.
      \item The \sQuote{MANUAL} help pages were updated.
      \item Confidence intervals (CIs) are renamed to reference intervals (RIs) in the \sQuote{MANUAL} and legends for being a more appropriate definition for these intervals (no argument is modified).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Functions \code{N.plot} and \code{p.plot} displayed all four types of records in the legend, even when only one type of record was plotted. Now the legend only shows the type of records that are plotted.
    }
  }
}

\section{Changes in Version 2.1.0 (05 August 2021)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function \code{change.point} to perform the Castillo-Mateo change-point detection tests was added. 
      \item Previously, if the series had \code{NA} values, all the record indicators from that time took \code{NA} value. This made the functions return \code{NA} p-values, etc. Now the \code{NA} values are treated as non-records, i.e., \eqn{I_t = 0}, unless they are the first observation, i.e., \eqn{I_1 = 1}.
      \item \code{records} function only allowed to show points, now it also shows lines.
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item The three missing values of the \code{TX_Zaragoza} dataset now take the value \code{NA}, in previous Versions these values were imputed. This change does not affect the \code{ZaragozaSeries} dataset.
      \item The \sQuote{MANUAL} help pages were updated.
      \item The \sQuote{VIGNETTES} file was updated.
    }
  }
}

\section{Changes in Version 2.0.0 (13 February 2021)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function \code{brown.method} to perform the the Brown's method in the (weighted)-number-of-records was added.
      \item New function \code{fisher.method} to perform Fisher's method was added.
      \item New function \code{R.record} to compute record values.
      \item New function \code{rcrm} to generate random samples of the classical record model.
      \item New function \code{series_ties} to summarize the presence of record ties (or weak records) mainly due to rounding effects.
      \item New dataset \code{Olympic_records_200m} with record times and record values of the Olympic 200-meter race.
      \item New function \code{series_record} to apply the tools within \pkg{RecordTest} if only the record times are available and not the complete series.
      \item New argument \code{weak} if weak records are also computed in functions that compute the record random variables: \code{I.record}, \code{L.record}, \code{N.record}, etc.
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item The \sQuote{MANUAL} help pages were updated.
      \item The \sQuote{VIGNETTES} file was updated.
      \item All functions have undergone a change in their arguments, providing new features in most cases.
      \item The functions \code{foster.test} and \code{foster.plot} now admit more statistics based on different combinations of records.
      \item The function \code{series_uncor} has change the criteria of considering the first and last columns uncorrelated, thus slightly different results can come from previous Versions. As new features it admits more types of uncorrelation and more arguments.
      \item The dataset \code{TX_Zaragoza} has been updated with more data.
      \item The matrix \code{ZaragozaSeries} has more rows and one more column according to the changes in \code{TX_Zaragoza} and \code{series_uncor}.
      \item Functions have been joined or renamed as follows:
      \itemize{
        \item From \code{L_global.test} renamed to \code{global.test}, and now a global statistic for different types of record can be computed for the tests with two-sided alternatives \code{p.chisq.test}, \code{p.regression.test}, \code{lr.test} and \code{score.test}.
        \item From \code{N_normal.test}, \code{N_exactPB.test} and \code{P_exactPB.test} joint to \code{N.test}, and new features added.
        \item From \code{P_regression.test} and \code{P_regression.plot} renamed to \code{p.regression.test} and \code{p.plot}, and new features added.
        \item From \code{P_chisq.test} renamed to \code{p.chisq.test}.
        \item From \code{L_lr.test} and \code{L_lm.test} renamed to \code{lr.test} and \code{score.test}.
        \item From \code{M.record} renamed to \code{S.record}.
      }
    }
  }  
  \subsection{BUG FIXES}{
    \itemize{
      \item The p-values in Monte-Carlo and Poisson-Binomial tests were computed with the expression \eqn{P(X > X_0)} instead of \eqn{P(X \ge X_0)} in all previous Versions, this problem was fixed.
      \item Function \code{foster.plot} showed the evolution of the \eqn{D} statistic improperly, this problem was fixed.
    }
  }
}

\section{Changes in Version 1.0.1 (05 October 2020)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The \sQuote{README} and \sQuote{VIGNETTES} files were added.
      \item New function \code{foster.test} to perform the Foster-Stuart and Diersen-Trenkler tests was added. 
      \item New function \code{foster.plot} to perform plots based on the Foster-Stuart and Diersen-Trenkler statistics was added. 
      \item New function \code{L_global.test} to perform likelihood ratio and Lagrange multiplier tests on record occurrence was added.
      \item New function \code{L_lm.test} to perform Lagrange multiplier test on record occurrence was added.
      \item New function \code{series_untie} to deal with record ties was added.
      \item Support for parallel computing was added in Monte Carlo functions: \code{L_global.test}, \code{L_lm.test}, \code{L_lr.test}.
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item Functions removed: \code{dailymean.fun}, \code{dailysd.fun}, \code{N_joint.plot}, \code{N_joint.test}, \code{std.fun}.
      \item Functions have been joined or renamed as follows:
      \itemize{
        \item From \code{I.rec} and \code{I.rec.matrix} joined to the S3 method \code{I.record}.
        \item From \code{L.rec} and \code{L.rec.matrix} joined to \code{L.record}.
        \item From \code{M.rec} and \code{P.rec} renamed to \code{M.record} and \code{P.record}, respectively.
        \item From \code{N.rec} and \code{N.rec.matrix} joined to \code{N.record}.
        \item From \code{Nmean.rec} renamed to \code{Nmean.record}.
        \item From \code{DaySeries.fun}, \code{double.fun}, \code{mintime.cor.test}, \code{value.rec} renamed to \code{series_split}, \code{series_double}, \code{series_uncor}, \code{records}, respectively.
        \item From \code{L_lr.test} and \code{L_MC.test} joined to \code{L_lr.test}.
      }
      \item Data \code{ZaragozaSeries}: now the \code{series_double} function has not been applied to it. 
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item The p-value and statistic of \code{P_regression.test} was calculated incorrectly in Version 0.1.2, this problem was fixed.
    }
  }
}

\section{Changes in Version 0.1.2 (09 February 2020)}{
  \subsection{OTHER CHANGES}{
    \itemize{
      \item Update by CRAN policies in some \sQuote{Imports} packages.
      \item Minor changes.
    }
  }
}