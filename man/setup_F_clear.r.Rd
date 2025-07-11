% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_clearance.R
\name{setup_F_clear.r}
\alias{setup_F_clear.r}
\title{Set up default F_clear}
\usage{
\method{setup_F_clear}{r}(mod_clear, r_opts, Xpars, nStrata)
}
\arguments{
\item{mod_clear}{is a string that defines the S3 class}

\item{r_opts}{a set of options to override the defaults}

\item{Xpars}{an \code{Xmod} object}

\item{nStrata}{the number of population strata}
}
\value{
an \code{Xmod} object
}
\description{
Implements \link{F_clear} for the SIP model.
}
