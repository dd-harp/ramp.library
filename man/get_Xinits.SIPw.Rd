% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-xde-SIPw.R
\name{get_Xinits.SIPw}
\alias{get_Xinits.SIPw}
\title{Return initial values as a vector from a SIPw model}
\usage{
\method{get_Xinits}{SIPw}(pars, i)
}
\arguments{
\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
none
}
\description{
This method dispatches on the type of \code{pars$Xpar[[i]]}
}
