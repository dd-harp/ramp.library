% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SIRS.R
\name{set_Xinits.SIRS}
\alias{set_Xinits.SIRS}
\title{Return the parameters as a list}
\usage{
\method{set_Xinits}{SIRS}(pars, i = 1, Xopts = list())
}
\arguments{
\item{pars}{an \strong{\code{xds}} object}

\item{i}{the vector species index}

\item{Xopts}{a named list}
}
\value{
an **`xds`** object
}
\description{
This method dispatches on the type of `pars$Xpar[[i]]`.
}
