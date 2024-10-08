% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SIRS.R
\name{get_Xinits.SIRS}
\alias{get_Xinits.SIRS}
\title{Return initial values as a vector}
\usage{
\method{get_Xinits}{SIRS}(pars, i)
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
