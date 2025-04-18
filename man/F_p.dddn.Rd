% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-bionomics-dddn.R
\name{F_p.dddn}
\alias{F_p.dddn}
\title{Dawn, day, dusk, night model for the human fraction}
\usage{
\method{F_p}{dddn}(t, vars, p_par)
}
\arguments{
\item{t}{current simulation time}

\item{vars}{exogenous variables}

\item{p_par}{a \link{list}}
}
\value{
a \link{numeric} vector of length \code{nPatches}
}
\description{
Implements \link{F_p} for a dddn model
}
