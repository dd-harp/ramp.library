% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/host-xde-SIR.R
\name{dXdt.SIR}
\alias{dXdt.SIR}
\title{Compute the derivatives for parasite infection dynamics in human population strata}
\usage{
\method{dXdt}{SIR}(t, y, pars, i)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a \link{numeric} vector
}
\description{
Implements \link{dXdt} for the SIR model
}
