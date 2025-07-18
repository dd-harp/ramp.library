% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria_detect.R
\name{F_detect_vars}
\alias{F_detect_vars}
\title{Set up a model for detection}
\usage{
F_detect_vars(vars, Xpars, mod_detect)
}
\arguments{
\item{vars}{variables as a named list}

\item{Xpars}{an \strong{X} model object}

\item{mod_detect}{an detection model object}
}
\value{
a \link{numeric} vector of length \code{nStrata}
}
\description{
A model could have models
for detection by different methods, so \code{F_detect}
must pass the specific detection model.
}
