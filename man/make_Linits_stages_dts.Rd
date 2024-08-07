% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L-dts-stages.R
\name{make_Linits_stages_dts}
\alias{make_Linits_stages_dts}
\title{Make inits for stages_dts competition aquatic mosquito model}
\usage{
make_Linits_stages_dts(
  nHabitats,
  Lopts = list(),
  L10 = 1,
  L20 = 0,
  L30 = 0,
  L40 = 0,
  P0 = 0
)
}
\arguments{
\item{nHabitats}{the number of habitats in the model}

\item{Lopts}{a [list] that overwrites default values}

\item{L10}{initial conditions}

\item{L20}{initial conditions}

\item{L30}{initial conditions}

\item{L40}{initial conditions}

\item{P0}{initial conditions}
}
\value{
a [list] with Linits added
}
\description{
Make inits for stages_dts competition aquatic mosquito model
}
