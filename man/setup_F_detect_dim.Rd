% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_detect.R
\name{setup_F_detect_dim}
\alias{setup_F_detect_dim}
\title{Setup an \code{F_detect.dim} model object}
\usage{
setup_F_detect_dim(c0, pda, nStrata)
}
\arguments{
\item{c0}{probability of detection per detectious cite, no immunity}

\item{pda}{a model for detection by aoi, motivated by a model for parasite densities and detection}

\item{nStrata}{the number of population strata}
}
\value{
an \code{F_detect.c} model object
}
\description{
Set up a F_detect model
with immunity of the form \deqn{d_0(1- e^{-\alpha}}
}
