#!/usr/bin/gnuplot
#  AMDG
# JMJ-V!
# FdS 24.01.2020
# 
# The function Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,warmPop) returns the fit,
# given in Aoyama, Marleau, Mordasini & Ikoma (2020),
# to the radius as a function of accretion rate Mdot and planet mass Mp
# for the cold- (Pop='cold') or warm-start (Pop='warm') population.
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
# To test the function, you can use the macros in the following way
# (ACHTUNG: will call "reset" first!):
# 
# gnuplot> load "./Rp_from_Mdot_Mp_Popsynthfit.gnu"
# gnuplot> set macros
# gnuplot> @test_Rp_fit_Mdot_warm
# gnuplot> @test_Rp_fit_Mp_cold
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
#  Note: use real() to be able to use the function at Mp < 1
#  for esthetic reasons
# 
Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,Pop) = ( lm2 = log10(Mdot/1e-2), \
        Pop eq 'warm' ? real(0.411-0.244*lm2+3.45*exp(0.762*lm2) + (-0.489-0.0961*lm2+0.652*exp(0.353*lm2))*(Mp-1.0) + (-0.228-0.00106*lm2+0.226*exp(0.000220*lm2))*(Mp-1.0)**2.) \
                     : real(1.53+0.111*lm2+1.06*exp(0.906*lm2) + (-0.195-0.0307*lm2+0.0977*exp(0.000695*lm2))*(Mp-1.0) + (-0.250+0.000276*lm2+0.254*exp(0.000214*lm2))*(Mp-1.0)**2.) \
                             )

# simple test functions
# 
test_Rp_fit_Mdot_warm = "reset; set title 'WARM POPULATION'; set key top left Left rever title 'Mass (M_J)'; set aut; set log x; set xlabel 'dM/dt (M_E/yr)'; set ylabel 'Radius (R_J)'; plot [1e-5:1e-2] for [i in '20 15 10 5 3 1'] Rp_from_Mdot_Mp_Popsynthfit(x,i,'warm') w l lw 2 t i"
test_Rp_fit_Mp_warm = "reset; set title 'WARM POPULATION'; set key top left Left rever title 'dM/dt (M_E/yr)'; set aut; set xlabel 'Mass (M_J)'; set ylabel 'Radius (R_J)'; plot [0:30] for [i in '1e-1 2e-2 1e-2 1e-3 1e-4 1e-5'] Rp_from_Mdot_Mp_Popsynthfit(i,x,'warm') w l lw 2 t i"

print "# Function Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,Pop) loaded"
