\name{NEWS}
\title{News for Package \pkg{RecordTest}}
\newcommand{\cpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\issue}{\href{https://github.com/JorgeCastilloMateo/RecordTest/issues/#1}{(issue #1)}}

\section{Changes to RecordTest Version 1.0.0 (14 July 2020)}{
  \itemize{
    \item Big changes and improvements. Most functions have undergone a change in their functionalities.
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
    \item New functions: \code{foster.plot}, \code{foster.test}, \code{L_global.test}, \code{L_lm.test}, \code{series_untie}.
    \item Data \code{ZaragozaSeries}: now the \code{series_double} function has not been applied to it.
    \item Problems solved: the p-value and statistic of \code{P_regression.test} was calculated incorrectly in version 0.1.2, this problem has been solved.
    \item Support for parallel computing added in Monte Carlo functions: \code{L_global.test}, \code{L_lm.test}, \code{L_lr.test}.
  }
}

\section{Changes to RecordTest Version 0.1.2 (9 February 2020)}{
  \itemize{
    \item Small changes. They do not affect the user experience.
    \item Updates for CRAN policies on Imports packages.
  }
}