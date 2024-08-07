% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MYZ-xde-RMG.R
\name{F_fqZ.RMG_xde}
\alias{F_fqZ.RMG_xde}
\title{Blood feeding rate of the infective mosquito population}
\usage{
\method{F_fqZ}{RMG_xde}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{a \link{list}}

\item{s}{the species index}
}
\value{
a [numeric] vector of length `nPatches`
}
\description{
Implements [F_fqZ] for the RMG_xde model.
}
