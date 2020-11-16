_Saint-Maurice, commander of the Theban legion_

# `St-Moritz`: Suite of Tools to Model Observations of accRetIng planeTZ

This collects a few simple functions or scripts that can be useful for studying accreting (forming) planets.

## Contents

1. `Rp_from_Mdot_Mp_Popsynthfit.*`: Planet radius as a function of accretion rate and mass, from the population syntheses of [Mordasini et al. (2012)](http://adsabs.harvard.edu/abs/2012A%26A...547A.112M), as fit by [Aoyama, Marleau, Mordasini & Ikoma (2020, subm.)](https://arxiv.org/abs/2011.06608). The files are organised by language (`gnuplot` (4.2/5), `idl` (8.7), `python` (2.7); other versions might/should work too).

1. `accretionflowproperties`: Computes approximate properties of the accretion flow onto a forming gas giant, including at the shock: density, temperature, velocity as a function of radius and time since the beginning of the fall. Based on the (semi-)analytical and simulation results of [Marleau et al. (2019b)](https://ui.adsabs.harvard.edu/abs/2019ApJ...881..144M). This submodule can be found under [https://github.com/gabrielastro/accretionflowproperties](https://github.com/gabrielastro/accretionflowproperties).

## Notes

Author: Gabriel-Dominique Marleau (c) 2020

Collaborator: [Yuhiko Aoyama](http://www.aoyama.saloon.jp/)

Contact: `uni-tuebingen.de` with `gabriel.marleau` and `@` in front

Comments and requests welcome!
