% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_transmit.R
\name{F_transmit_y}
\alias{F_transmit_y}
\title{Compute Infectiousness from the State Variables}
\usage{
F_transmit_y(y, Xpars)
}
\arguments{
\item{y}{state vector}

\item{Xpars}{an \code{Xmod} object}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
This method computes infectiousness,
the probability a mosquito would become infected
after blood feeding on a human, by strata
}
