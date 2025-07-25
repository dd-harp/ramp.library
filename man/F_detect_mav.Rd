% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_detect.R
\name{F_detect_mav}
\alias{F_detect_mav}
\title{Infection blocking pre-erythrocytic immunity}
\usage{
F_detect_mav(vh, moi, aoi, d0, imm_rt, pd_aoi)
}
\arguments{
\item{vh}{cumulative exposure}

\item{moi}{mean multiplicity of infection}

\item{aoi}{mean age of infection}

\item{d0}{maximum probability of detection}

\item{imm_rt}{effect modivication by}

\item{pd_aoi}{parasite densities by aoi}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
Implements \link{F_detect_vars.mav} for the SIP model.
}
