% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SEISd.R
\name{get_Xinits.SEISd}
\alias{get_Xinits.SEISd}
\title{Return initial values as a vector}
\usage{
\method{get_Xinits}{SEISd}(pars, i = 1)
}
\arguments{
\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a \link{numeric} vector
}
\description{
This method dispatches on the type of \code{pars$Xpar}.
}
