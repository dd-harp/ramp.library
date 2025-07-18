% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_transmit.R
\name{setup_F_transmit}
\alias{setup_F_transmit}
\title{Setup Infectiousness}
\usage{
setup_F_transmit(mod_c, c_opts, Xpars, nStrata)
}
\arguments{
\item{mod_c}{is a string that defines the S3 class}

\item{c_opts}{a set of options to override the defaults}

\item{Xpars}{an \code{Xmod} object}

\item{nStrata}{the number of population strata}
}
\value{
an \code{Xmod} model object
}
\description{
Set up an \code{F_transmit} model object
}
