% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-xde-RMG.R
\name{make_MYZpar_RMG}
\alias{make_MYZpar_RMG}
\title{Make parameters for RM ODE adult mosquito model}
\usage{
make_MYZpar_RMG(
  nPatches,
  MYZopts = list(),
  eip = 12,
  g = 1/12,
  sigma_b = 1/8,
  sigma_q = 1/8,
  mu = 0,
  f = 0.5,
  q = 0.95,
  nu = 1,
  eggsPerBatch = 60
)
}
\arguments{
\item{nPatches}{is the number of patches, an integer}

\item{MYZopts}{a \link{list} of values that overwrites the defaults}

\item{eip}{extrinsic incubation period}

\item{g}{mosquito mortality rate}

\item{sigma_b}{emigration rate while blood feeding}

\item{sigma_q}{emigration rate while egg laying}

\item{mu}{emigration loss}

\item{f}{feeding rate}

\item{q}{human blood fraction}

\item{nu}{oviposition rate, per mosquito}

\item{eggsPerBatch}{eggs laid per oviposition}
}
\value{
a \link{list}
}
\description{
Make parameters for RM ODE adult mosquito model
}
