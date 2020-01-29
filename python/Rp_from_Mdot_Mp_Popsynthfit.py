#  AMDG
# JMJ-V!
# FdS 24.01.2020
# 
# The function Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,warmPop) returns the fit,
# given in Aoyama, Marleau, Mordasini & Ikoma (2020),
# to the radius as a function of accretion rate Mdot and planet mass Mp
# for the cold- (Pop='cold') or warm-start (Pop='warm') populations
# of Mordasini et al. (2012).
# The fit is valid for approximately (please see the paper):
# 
#          1 MJ <  Mp  < 30 MJ
#    1e-5 ME/yr < Mdot < 2e-2 ME/yr
# 
# -- Input --
#   Mdot: in ME/yr = 1.893e20 g/s
#     Mp: in MJ = 1.898e30 g
# 
# -- Output --
#     RP: in RJ = 7.15e9 cm
# 
# If you use this function, please cite: Aoyama et al. (2020). Thank you!
# 
# Comments and requests welcome!
# 
# Gabriel-Dominique Marleau (c) February 2020
# gabriel.marleau@uni-tuebingen.de, gabriel.marleau@space.unibe.ch
# 
# -------------------------------------------------------------
# 
import numpy as np

def Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,Pop='warm'):
	lm2 = np.log10(Mdot/1e-2)
	if Pop == 'warm':
		return 0.411-0.244*lm2+3.45*np.exp(0.762*lm2) \
		      + (-0.489-0.0961*lm2+0.652*np.exp(0.353*lm2))*(Mp-1.0) \
		      + (-0.228-0.00106*lm2+0.226*np.exp(0.000220*lm2))*(Mp-1.0)**2.
	else:
		return 1.53+0.111*lm2+1.06*np.exp(0.906*lm2) \
		      + (-0.195-0.0307*lm2+0.0977*np.exp(0.000695*lm2))*(Mp-1.0) \
		      + (-0.250+0.000276*lm2+0.254*np.exp(0.000214*lm2))*(Mp-1.0)**2.
