% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L-dts-stages.R
\name{dLdt.stages_dts}
\alias{dLdt.stages_dts}
\title{Derivatives for aquatic stage mosquitoes}
\usage{
\method{dLdt}{stages_dts}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \strong{\code{xds}} object}

\item{s}{the species index}
}
\value{
a [numeric] vector
}
\description{
Implements [dLdt] for the stages_dts competition model.
}
