% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_infect.R
\name{F_infect_y.im}
\alias{F_infect_y.im}
\title{Infections per Infectious Bite}
\usage{
\method{F_infect_y}{im}(y, Xpars)
}
\arguments{
\item{y}{state vector}

\item{Xpars}{an \code{Xmod} object}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
The basic model for the
risk of infection, per infectious bite:
a constant probability of
infection, per infectious bite
}
