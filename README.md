# ramp.library <br><br> A Model Library for RAMP

This package -- `ramp.library` -- contains reusable code that has been rigorously tested and that implements a large number of dynamical model families and other algorithms taken from the literature describing malaria and other mosquito-transmitted pathogens (see Reiner, *et al.* 2013)^[Reiner RC Jr, Perkins TA, Barker CM, Niu T, Chaves LF, Ellis AM, et al. A systematic review of mathematical models of mosquito-borne pathogen transmission: 1970-2010. J R Soc Interface. 2013;10: 20120921.]. The supporting code was designed to be modular, and plug-and-play. The modular design makes it possible to break down published models to serve as the dynamical components in new models for malaria. The model families in `ramp.library` supports nimble model building that can be used by:

+ [`ramp.xde`](https://github.com/dd-harp/ramp.xde) implements a modular, flexible, and extensible framework for building systems of ordinary and delay differential equations models for malaria and other mosquito-transmitted pathogens. The mathematic framework for `ramp.xde` was explained in [Spatial Dynamics of Malaria Transmission](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010684){target="_blank"}
^[Wu SL, Henry JM, Citron DT, Mbabazi Ssebuliba D, Nakakawa Nsumba J, Sánchez C. HM, et al. (2023) Spatial dynamics of malaria transmission. PLoS Comput Biol 19(6): e1010684. https://doi.org/10.1371/journal.pcbi.1010684] 

+ [`ramp.dts`](https://github.com/dd-harp/ramp.dts) implements a modular, flexible, and extensible framework for building discrete time systems, including both deterministic and stochastic difference equations for malaria and other mosquito-transmitted pathogens.  
