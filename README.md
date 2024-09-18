# **`ramp.library`** - A Code Library for [**`ramp.xds`**](https://dd-harp.github.io/ramp.xds/) 
## *Dynamical Systems for the Epidemiology, Spatial Transmission Dynamics, and Control of Malaria and Other Mosquito-borne Pathogens*

<!-- badges: start -->
[![PLoS Computational Biology](https://img.shields.io/badge/doi%3A10.1371%2Fjournal.pcbi.1010684-cornflowerblue?style=for-the-badge&logoColor=cornflowerblue&label=PLoS%20Comp%20Bio&labelColor=slategrey&link=https%3A%2F%2Fjournals.plos.org%2Fploscompbiol%2Farticle%3Fid%3D10.1371%2Fjournal.pcbi.1010684)](https://doi.org/10.1371/journal.pcbi.1010684)

<!-- badges: end -->


## Installation 

**`ramp.library`** is a code library that extends [**`ramp.xds`.**](https://dd-harp.github.io/ramp.xds/) It aims to hold a comprehensive set of published models with code that has been verified and tested. 

To install the latest version of **`ramp.xds`** from GitHub, run the following lines of code in an R session.

```
library(devtools)
devtools::install_github("dd-harp/ramp.xds")
```

To install the latest version of **`ramp.library`** from Github, run the following line in an R session: 
```
devtools::install_github("dd-harp/ramp.library")
```

## What is RAMP?

RAMP -- **R**obust **A**nalytics for **M**alaria **P**olicy -- is a bespoke inferential system for malaria decision support and adaptive malaria control. A core goal for RAMP is to characterize, quantify, and propagate uncertainty in conventional analysis and through simulation-based analytics.

## What is **`ramp.xds`**?

**`ramp.xds`** is an R software package that supports nimble model building for simulation-based analytics and malaria research. It was designed to help research scientists and policy analysts set up, analyze, solve, and apply dynamical systems models describing the epidemiology, spatial transmission dynamics, and control of malaria and other mosquito-transmitted pathogens. The software also supports nimble model building and analysis for mosquito ecology, with the capability to handle forcing by weather and other exogenous variables. 

The software was designed around a rigorous mathematical framework for modular model building, described in [Spatial Dynamics of Malaria Transmission](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010684) (Wu SL, *et al.* 2023. PLoS Computational Biology)^[Wu SL, Henry JM, Citron DT, Mbabazi Ssebuliba D, Nakakawa Nsumba J, Sánchez C. HM, et al. (2023) Spatial dynamics of malaria transmission. PLoS Comput Biol 19(6): e1010684. https://doi.org/10.1371/journal.pcbi.1010684]. The mathematical framework has now been extended to cover *exogenous forcing* by weather and vector control. 

## What is **`ramp.library`**?

**`ramp.library`** is a code library with models that can be used by **`ramp.xds`.** Thee two packages are part of a suite of R packages developed to support RAMP: 

+ [**`ramp.xds`.**](https://dd-harp.github.io/ramp.xds/) is the core computational engine for simulation-based analytics. It includes a basic set of models -- enough to design, verify, and demonstrate the basic features of modular software. 

+  **`ramp.library`** is an extended library of stable code that has been tested and verified. It includes a large set of model families published in peer review that are not included in **`ramp.xds`** The ability to reuse code reduces the costs of replicating studies. Through this library, **`ramp.xds`** also supports nimble model building and analytics for other mosquito-borne pathogens. 

+ [**`ramp.work`**](https://dd-harp.github.io/ramp.work/) includes algorithms to apply the framework, include code to fit models to data and to do constrained optimization 

+ [**`ramp.malaria`**](https://dd-harp.github.io/ramp.malaria/) includes a large set of models illustrating capabilities of **`ramp.xds`** 

**`ramp.library`**  is under active development. 

# **`ramp.library`** <br><br> A Model Library for RAMP

This package -- `ramp.library` -- contains reusable code that has been rigorously tested and that implements a large number of dynamical model families and other algorithms taken from the literature describing malaria and other mosquito-transmitted pathogens (see Reiner, *et al.* 2013)^[Reiner RC Jr, Perkins TA, Barker CM, Niu T, Chaves LF, Ellis AM, et al. A systematic review of mathematical models of mosquito-borne pathogen transmission: 1970-2010. J R Soc Interface. 2013;10: 20120921.]. The supporting code was designed to be modular, and plug-and-play. The modular design makes it possible to break down published models to serve as the dynamical components in new models for malaria. 


