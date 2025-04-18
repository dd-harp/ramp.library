% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-xde-SIPw.R
\name{parse_Xorbits.SIPw}
\alias{parse_Xorbits.SIPw}
\title{Parse the output of deSolve and return variables for SIPw models}
\usage{
\method{parse_Xorbits}{SIPw}(outputs, pars, i)
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
Implements \link{parse_Xorbits} for SIPw models
}
