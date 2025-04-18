% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/human-SIP.R
\name{F_ni.SIP}
\alias{F_ni.SIP}
\title{Compute the "true" prevalence of infection / parasite rate}
\usage{
\method{F_ni}{SIP}(vars, Xpar)
}
\arguments{
\item{vars}{a list with the variables attached by name}

\item{Xpar}{a list defining a model for human}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_ni} for the SIP model.
}
