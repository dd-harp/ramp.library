% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-xde-SIPw.R
\name{F_b.SIPw}
\alias{F_b.SIPw}
\title{Infection blocking pre-erythrocytic immunity}
\usage{
\method{F_b}{SIPw}(y, pars, i)
}
\arguments{
\item{y}{state vector}

\item{pars}{an \strong{\code{xds}} object}

\item{i}{the host species index}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_b} for SIPw models
}
