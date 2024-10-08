% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MYZ-xde-RMG.R
\name{F_eggs.RMG}
\alias{F_eggs.RMG}
\title{Number of eggs laid by adult mosquitoes}
\usage{
\method{F_eggs}{RMG}(t, y, pars, s)
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
Implements [F_eggs] for the RMG model.
}
