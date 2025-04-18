% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-dts-ENBRQ.R
\name{F_eggs.ENBRQ_dts}
\alias{F_eggs.ENBRQ_dts}
\title{Number of eggs laid by adult mosquitoes}
\usage{
\method{F_eggs}{ENBRQ_dts}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{a \link{list}}

\item{s}{the species index}
}
\value{
a \link{numeric} vector of length \code{nPatches}
}
\description{
Implements \link{F_eggs} for the ENBRQ_dts model.
}
