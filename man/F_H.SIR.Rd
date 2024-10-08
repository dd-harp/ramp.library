% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SIR.R
\name{F_H.SIR}
\alias{F_H.SIR}
\title{Size of effective infectious human population}
\usage{
\method{F_H}{SIR}(t, y, pars, i)
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
Implements [F_H] for the SIR model.
}
