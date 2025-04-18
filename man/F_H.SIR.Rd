% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/host-xde-SIR.R
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
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_H} for the SIR model.
}
