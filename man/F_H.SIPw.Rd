% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-SIPw.R
\name{F_H.SIPw}
\alias{F_H.SIPw}
\title{Size of effective infectious human population}
\usage{
\method{F_H}{SIPw}(t, y, pars, i)
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
Implements [F_X] for SIPw models
}
