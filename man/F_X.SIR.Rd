% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SIR.R
\name{F_X.SIR}
\alias{F_X.SIR}
\title{Size of effective infectious human population}
\usage{
\method{F_X}{SIR}(t, y, pars, i)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a [numeric] vector of length `nStrata`
}
\description{
Implements [F_X] for the SIS model.
}
