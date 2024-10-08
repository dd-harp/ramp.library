% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/X-xde-SEIRV.R
\name{make_Xinits_SEIRV}
\alias{make_Xinits_SEIRV}
\title{Make initial values for the SEIRV human model, with defaults}
\usage{
make_Xinits_SEIRV(nStrata, H, Xopts = list(), I = 1, E = 0, R = 1, V = 1)
}
\arguments{
\item{nStrata}{the number of strata in the model}

\item{H}{the initial value for H}

\item{Xopts}{a [list] to overwrite defaults}

\item{I}{the initial value for I}

\item{E}{the initial value for E}

\item{R}{the initial values for R}

\item{V}{the initial values for V}
}
\value{
a [list]
}
\description{
Make initial values for the SEIRV human model, with defaults
}
