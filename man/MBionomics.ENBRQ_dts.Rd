% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-dts-ENBRQ.R
\name{MBionomics.ENBRQ_dts}
\alias{MBionomics.ENBRQ_dts}
\title{Reset bloodfeeding and mortality rates to baseline}
\usage{
\method{MBionomics}{ENBRQ_dts}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \code{xds} object}

\item{s}{the species index}
}
\value{
the model as a \link{list}
}
\description{
Implements \link{MBionomics} for the ENBRQ_dts model
}
