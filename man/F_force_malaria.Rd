% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-forcing.R
\name{F_force_malaria}
\alias{F_force_malaria}
\title{Mass Treatment}
\usage{
F_force_malaria(t, y, Xpars, forcing_mod)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{Xpars}{an \code{Xmod} object}

\item{forcing_mod}{a function}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Set up a function to force some aspect of malaria
}
