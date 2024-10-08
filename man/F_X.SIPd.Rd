% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SIPd.R
\name{F_X.SIPd}
\alias{F_X.SIPd}
\title{Size of effective infectious human population}
\usage{
\method{F_X}{SIPd}(t, y, pars, i)
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
Implements [F_X] for the SIPd model.
}
