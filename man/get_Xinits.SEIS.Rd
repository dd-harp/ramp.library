% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SEIS.R
\name{get_Xinits.SEIS}
\alias{get_Xinits.SEIS}
\title{Return initial values as a vector}
\usage{
\method{get_Xinits}{SEIS}(pars, i)
}
\arguments{
\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a [numeric] vector
}
\description{
This method dispatches on the type of `pars$Xpar`.
}
