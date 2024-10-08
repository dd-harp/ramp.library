% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SIR.R
\name{F_b.SIR}
\alias{F_b.SIR}
\title{Infection blocking pre-erythrocytic immunity}
\usage{
\method{F_b}{SIR}(y, pars, i)
}
\arguments{
\item{y}{state vector}

\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a [numeric] vector of length `nStrata`
}
\description{
Implements [F_b] for the SIR model.
}
