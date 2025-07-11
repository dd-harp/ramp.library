% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/host-xde-SEIR.R
\name{F_prevalence.SEIR}
\alias{F_prevalence.SEIR}
\title{Compute the "true" prevalence of infection / parasite rate}
\usage{
\method{F_prevalence}{SEIR}(vars, Xpar)
}
\arguments{
\item{vars}{a list with the variables attached by name}

\item{Xpar}{a list defining a model for human}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_prevalence} for the SEIR model.
}
