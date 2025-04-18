% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SIP.R
\name{add_lines_X_SIP}
\alias{add_lines_X_SIP}
\title{Add lines for the density of infected individuals for the SIP model}
\usage{
add_lines_X_SIP(
  time,
  XH,
  nStrata,
  clrs = c("darkblue", "darkred", "darkgreen"),
  llty = 1
)
}
\arguments{
\item{time}{time points for the observations}

\item{XH}{a list with the outputs of parse_Xorbits.SIP}

\item{nStrata}{the number of population strata}

\item{clrs}{a vector of colors}

\item{llty}{an integer (or integers) to set the \code{lty} for plotting}
}
\description{
Add lines for the density of infected individuals for the SIP model
}
