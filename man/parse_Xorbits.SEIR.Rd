% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/host-xde-SEIR.R
\name{parse_Xorbits.SEIR}
\alias{parse_Xorbits.SEIR}
\title{Parse the output of deSolve and return variables for the SEIR model}
\usage{
\method{parse_Xorbits}{SEIR}(outputs, pars, i)
}
\arguments{
\item{outputs}{a \link{matrix} of outputs from deSolve}

\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
none
}
\description{
Implements \link{parse_Xorbits} for the SEIR model
}
