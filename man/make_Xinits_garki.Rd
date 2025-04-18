% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-garki.R
\name{make_Xinits_garki}
\alias{make_Xinits_garki}
\title{Make inits for garki human model.}
\usage{
make_Xinits_garki(
  nStrata,
  H,
  Xopts = list(),
  x2 = 0,
  y1 = 0,
  y2 = 0,
  y3 = 0,
  x3 = 0,
  x4 = 0
)
}
\arguments{
\item{nStrata}{is the number of population strata}

\item{H}{a \link{numeric} initial value for total human population density}

\item{Xopts}{a \link{list} with values to override default values}

\item{x2}{a \link{numeric} initial value for the variable x2}

\item{y1}{a \link{numeric} initial value for the variable y1}

\item{y2}{a \link{numeric} initial value for the variable y2}

\item{y3}{a \link{numeric} initial value for the variable y3}

\item{x3}{a \link{numeric} initial value for the variable x3}

\item{x4}{a \link{numeric} initial value for the variable x4}
}
\value{
none
}
\description{
Make inits for garki human model.
}
\note{
We use H instead of x, but other variables are passed in the same order as they are presented in the original paper.
}
