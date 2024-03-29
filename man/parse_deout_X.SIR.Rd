% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Xmod-SIR.R
\name{parse_deout_X.SIR}
\alias{parse_deout_X.SIR}
\title{Parse the output of deSolve and return variables for the SIR model}
\usage{
\method{parse_deout_X}{SIR}(deout, pars, i)
}
\arguments{
\item{deout}{a \link{matrix} of outputs from deSolve}

\item{pars}{a \link{list} that defines a model}

\item{i}{the host species index}
}
\value{
none
}
\description{
Implements [parse_deout_X] for the SIR model
}
