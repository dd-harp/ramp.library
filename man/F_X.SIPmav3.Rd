% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-SIPmav3.R
\name{F_X.SIPmav3}
\alias{F_X.SIPmav3}
\title{Compute Infectious Density, \eqn{X}, for \code{SIPmav3}}
\usage{
\method{F_X}{SIPmav3}(t, y, pars, i)
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
Implements \link{F_X} for SIPmav3 models
}
