% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_infect.R
\name{F_infect_vars.b}
\alias{F_infect_vars.b}
\title{Infections per Infectious Bite}
\usage{
\method{F_infect_vars}{b}(vars, Xpars)
}
\arguments{
\item{vars}{variables as a named list}

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
