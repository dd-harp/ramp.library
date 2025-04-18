% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-xde-SBRQ.R
\name{parse_MYZorbits.SBRQ}
\alias{parse_MYZorbits.SBRQ}
\title{Parse the output of deSolve and return variables for the SBRQ model}
\usage{
\method{parse_MYZorbits}{SBRQ}(outputs, pars, s)
}
\arguments{
\item{outputs}{a \link{matrix} of outputs from deSolve}

\item{pars}{a \link{list} that defines a model}

\item{s}{the species index}
}
\value{
none
}
\description{
Implements \link{parse_MYZorbits} for the SBRQ model
}
