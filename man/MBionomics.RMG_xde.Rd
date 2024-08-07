% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MYZ-xde-RMG.R
\name{MBionomics.RMG_xde}
\alias{MBionomics.RMG_xde}
\title{Reset bloodfeeding and mortality rates to baseline}
\usage{
\method{MBionomics}{RMG_xde}(t, y, pars, s)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{a \link{list}}

\item{s}{the species index}
}
\value{
a named [list]
}
\description{
Implements [MBionomics] for the RMG_xde model
}
