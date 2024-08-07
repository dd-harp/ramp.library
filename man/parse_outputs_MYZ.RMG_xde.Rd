% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MYZ-xde-RMG.R
\name{parse_outputs_MYZ.RMG_xde}
\alias{parse_outputs_MYZ.RMG_xde}
\title{Parse the output of deSolve and return variables for the RMG_xde model}
\usage{
\method{parse_outputs_MYZ}{RMG_xde}(outputs, pars, s)
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
Implements [parse_outputs_MYZ] for the RMG_xde model
}
