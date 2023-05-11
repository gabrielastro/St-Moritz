_Saint-Maurice, commander of the Theban legion_

# `St-Moritz`: Suite of Tools to Model Observations of accRetIng planeTZ

This collects a few simple functions or scripts that can be useful for studying accreting (forming) planets.

## Contents

1. `RpTeff_from_Mdot_Mp_Popsynthfit.gnu`, `Rp_from_Mdot_Mp_Popsynthfit.{pro,py}`: Planet radius and Teff as a function of accretion rate and mass, from the population syntheses of [Mordasini et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...547A.112M), as fit by [Aoyama, Marleau, Mordasini & Ikoma (2020)](https://arxiv.org/abs/2011.06608). The files are organised by language (`gnuplot` (4.2/5), `idl` (8.7), `python` (3.7); other versions might/should work too).  
**New (13.02.2021)**: The Teff fit and test functions have been added to all three languages.

1. `accretionflowproperties`: Computes approximate properties of the accretion flow onto a forming gas giant, including, at the shock, the density, temperature, velocity as a function of radius (distance from the planet) and time since the beginning of the fall. Assumes spherical symmetry. Based on the (semi-)analytical and simulation results of [Marleau et al. (2019b)](https://ui.adsabs.harvard.edu/abs/2019ApJ...881..144M). This submodule can be found under [https://github.com/gabrielastro/accretionflowproperties](https://github.com/gabrielastro/accretionflowproperties).

1. Upcoming, when time permits: the relations between accretion luminosity Lacc and line luminosity Lline from our [Marleau & Aoyama (2022) Research Note "Planetary Line-to-accretion Luminosity Scaling Relations: Extrapolating to Higher-order Hydrogen Lines"](https://ui.adsabs.harvard.edu/abs/2022RNAAS...6..262M), complementing the Lacc–Lline scalings in [Aoyama et al. (2021, ApJL)](http://ads.nao.ac.jp/abs/2021ApJ...917L..30A). The explicit relations given there are already available in `python` in Tomas Stolker's [`species` toolkit](https://github.com/tomasstolker/species); see the [tutorial in `species`](https://species.readthedocs.io/en/latest/tutorials/emission_line.html?highlight=accretion#Emission-line-analysis) analysing the Paschen beta emission from GQ Lup B as an example.

## Notes

Author: Gabriel-Dominique Marleau (c) 2023

Collaborator: [Yuhiko Aoyama](http://www.aoyama.saloon.jp/)

Contact: `uni-tuebingen.de` with `gabriel.marleau` and `@` in front

Note that the free, open-source [Gnu Data Language (`GDL`)](https://github.com/gnudatalanguage) and also the closed-sourced but free [Fawlty Language `FL`](https://www.flxpert.hu/fl/) exist as drop-in replacements for IDL!

Comments and requests welcome.
