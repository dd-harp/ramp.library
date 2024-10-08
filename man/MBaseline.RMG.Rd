% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MYZ-xde-RMG.R
\name{MBaseline.RMG}
\alias{MBaseline.RMG}
\title{Reset bloodfeeding and mortality rates to baseline}
\usage{
\method{MBaseline}{RMG}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \code{xds} object}

\item{s}{the species index}
}
\value{
a named [list]
}
\description{
Implements [MBaseline] for the RMG model
}
