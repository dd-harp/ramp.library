% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/malaria-SIPmav3.R
\name{make_Xpar_SIPmav3}
\alias{make_Xpar_SIPmav3}
\title{Make parameters for SIPmav3_xde human model, with defaults}
\usage{
make_Xpar_SIPmav3(
  nStrata,
  Xopts = list(),
  mod_b = "b",
  b_opts = list(b = 0.55),
  mod_c = "c",
  c_opts = list(c = 0.15),
  mod_r = "r",
  r_opts = list(r = 1/180),
  mod_lm = "d",
  lm_opts = list(d = 0.7),
  mod_pcr = "d",
  pcr_opts = list(d = 0.9),
  mod_rdt = "d",
  rdt_opts = list(d = 0.7),
  mod_severe_h = "p",
  severe_h_opts = list(p = 1/500),
  mod_severe_x = "rate",
  severe_x_opts = list(rate = 1/365/5),
  mod_moderate_h = "p",
  moderate_h_opts = list(p = 1/50),
  mod_moderate_x = "rate",
  moderate_x_opts = list(rate = 1/730),
  mod_mild_h = "p",
  mild_h_opts = list(p = 1/10),
  mod_mild_x = "rate",
  mild_x_opts = list(rate = 1/50),
  mod_treat = "p3",
  treat_opts = list(p_sev = 0.95, p_mod = 0.1, p_mild = 0.01),
  mod_F_xi = "xi",
  xi_opts = list(xi = 1/730),
  mass_treat_opts = list(),
  mass_test_treat_opts = list(),
  mod_seek_care = "xi",
  care_seeking_opts = list(xi = 1/730),
  mod_other_disease = "xi",
  other_disease_opts = list(xi = 1/730),
  eta = 1/25
)
}
\arguments{
\item{nStrata}{the number of population strata}

\item{Xopts}{a \link{list} that could overwrite defaults}

\item{mod_b}{a string, the S3 class name for F_infect}

\item{b_opts}{a list: options for setup_F_infect}

\item{mod_c}{a string, the S3 class name for F_transmit}

\item{c_opts}{a list: options for setup_F_transmit}

\item{mod_lm}{a string, the S3 class name for F_detect for light microscopy}

\item{lm_opts}{a list: options for setup_F_detect}

\item{mod_pcr}{a string, the S3 class name for F_detect for PCR}

\item{pcr_opts}{a list: options for setup_F_detect}

\item{mod_rdt}{a string, the S3 class name for F_detect for RDT}

\item{rdt_opts}{a list: options for setup_F_detect}

\item{severe_h_opts}{a list: options for setup_F_incidence_h}

\item{mod_severe_x}{a string, the S3 class name for F_incidence_x for severe malaria}

\item{severe_x_opts}{a list: options for setup_F_incidence_x}

\item{mod_moderate_h}{a string, the S3 class name for F_incidence_h for moderate malaria}

\item{moderate_h_opts}{a list: options for setup_F_incidence_h}

\item{mod_moderate_x}{a string, the S3 class name for F_incidence_x for moderate malaria}

\item{moderate_x_opts}{a list: options for setup_F_incidence_x}

\item{mod_mild_h}{a string, the S3 class name for F_incidence_h for mild malaria}

\item{mild_h_opts}{a list: options for setup_F_incidence_h}

\item{mod_mild_x}{a string, the S3 class name for F_incidence_x for mild malaria}

\item{mild_x_opts}{a list: options for setup_F_incidence_x}

\item{eta}{rate of loss of chemo-protection}

\item{mod_mod_severe_h}{a string, the S3 class name for F_incidence_h for severe malaria}

\item{c}{transmission probability (efficiency) from human to mosquito}

\item{r}{recovery rate}

\item{rho}{probability of successful treatment upon infection}

\item{sigma}{probability of treatment for treated, above background}

\item{xi}{background treatment rate}
}
\value{
a \link{list}
}
\description{
Make parameters for SIPmav3_xde human model, with defaults
}
