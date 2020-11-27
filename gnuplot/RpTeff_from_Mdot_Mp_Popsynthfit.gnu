#!/usr/bin/env gnuplot
#  AMDG
# JMJ-V!
# FdS 24.01.2020
# 
# This file contains two functions giving the radius and the Teff of forming planets.
# 
# A.)
# The function Rp_from_Mdot_Mp(Mdot,Mp,Pop) returns the fit,
# given in Aoyama et al. (2020),
# to the radius as a function of accretion rate Mdot and planet mass Mp
# for the cold- (Pop='cold') or warm-start (Pop='warm') populations
# of Mordasini et al. (2012).
# 
# B.)
# The function Teff_from_Mdot_Mp(Mdot,Mp,Tint,Pop,ffill) returns
# the fit of Aoyama et al. (2020, their §2.2.2) to the Teff
# of the accreting region of the planet's surface.
# It ensures approximate energy conservation. The details of the fit
# might change in future iterations of the code.
# 
# The non-accreting region is in principle emitting at Teff=Tint
# but this is not relevant for this code.
# 
# The fits are valid for approximately (please see the paper):
# 
#          1 MJ <  Mp  < 30 MJ
#    1e-5 ME/yr < Mdot < 2e-2 ME/yr
# 
# -- Input --
#   Mdot: in ME/yr = 1.893e20 g/s
#     Mp: in MJ = 1.898e30 g
# 
# -- Output --
#     Rp: in RJ = 7.15e9 cm
#   Teff: in K
# 
# Works with gnuplot >= 4.2* (or earlier)
# 
# To test the functions, you can use the macros in the following way
# (SIDE EFFECT: the test will call "reset"):
# 
# gnuplot> # silent = 1  ## uncomment to suppress messages after loading
# gnuplot> load "./Rp_from_Mdot_Mp.gnu"
# gnuplot> set macros
# gnuplot> @test_Rp_fit_Mdot_warm
# gnuplot> @test_Rp_fit_Mp_warm
# gnuplot> @test_Teff_fit_Mdot_warm
# gnuplot> @test_Teff_fit_Mp_warm
# 
# Each will produce a plot in the current terminal.
# 
# If you use any of this, please cite:
# 
#   Aoyama, Marleau, Mordasini & Ikoma (2020)
#   Spectral appearance of the planetary surface accretion shock:
#     Global spectra and hydrogen line profiles and fluxes
#   https://arxiv.org/abs/2011.06608
# 
# Thank you. Comments and requests welcome!
# 
# Gabriel-Dominique Marleau (c) 2020
# @uni-tuebingen.de with gabriel.marleau in front
# 
# -------------------------------------------------------------
# 

# print some information to screen after having loaded?
#   By default just three lines (verbose=1).
#   Suppress by setting "verbose=0" before "load"ing this file
if(!exists("verbose")) verbose = 1

#  Note: use real() to be able to use the function at Mp < 1
#  for esthetic reasons
# 
Rp_from_Mdot_Mp(Mdot,Mp,Pop) = ( lm2 = log10(Mdot/1e-2), \
        Pop eq 'warm' ? real(0.411-0.244*lm2+3.45*exp(0.762*lm2) + (-0.489-0.0961*lm2+0.652*exp(0.353*lm2))*(Mp-1.0) + (-0.228-0.00106*lm2+0.226*exp(0.000220*lm2))*(Mp-1.0)**2.) \
                     : real(1.53+0.111*lm2+1.06*exp(0.906*lm2) + (-0.195-0.0307*lm2+0.0977*exp(0.000695*lm2))*(Mp-1.0) + (-0.250+0.000276*lm2+0.254*exp(0.000214*lm2))*(Mp-1.0)**2.) \
                             )

Teff_from_Mdot_Mp(Mdot,Mp,Tint,Pop,ffill) = (Tint**4. + \
             fdown_Popsynthfit(log10(n_MdotMp(Mdot,Mp,Pop,ffill)),v_MdotMp(Mdot,Mp,Pop)/-1e5) \
             * Tacc_MdotMp(Mdot,Mp,Pop,ffill)**4.)**0.25

# -------------------------------------------------------------

# 
# the "accretion temperature" of the accreting region
#   - Here it is correct to divide by sigma_SB=5.67e-5 cgs and not 4*sigma_SB = a*c
#     because we want the effective temperatures, not the gas temperature
#     (see Marleau et al. 2019; http://adsabs.harvard.edu/abs/2019ApJ...881..144M)
#   - The rest of the planet's surface is emitting only with Tint
# 
Tacc_MdotMp(Mdot,Mp,Pop,ffill) = ( 6.67e-8*Mp*1.898e30*Mdot*1.893e20/ (ffill * 4*pi*(Rp_from_Mdot_Mp(Mdot,Mp,Pop)*7.15e9)**3. * 5.67e-5) )**0.25


# 
# fit to the data of Aoyama et al. (2018) and Aoyama et al. (2020; https://arxiv.org/abs/2011.06608)
# 
fdown1_Popsynthfit(lgn0,v0) = 0.703752-0.0967987*(lgn0-12)-0.0254579*abs(lgn0-12)**2.  +  (-0.00527886-0.00146833*(lgn0-12)-0.000321504*abs(lgn0-12)**2.)*(v0-100)  -9.91492e-06*abs(v0-100)**2.

# limit fdown to [0,1]
# 
## avoid possibly overwriting a user's function
StMoritz_min(a,b) = (a<b?a:b)
StMoritz_max(a,b) = (a>b?a:b)
fdown_Popsynthfit(lgn0,v0) = StMoritz_max( StMoritz_min(fdown1_Popsynthfit(lgn0,v0),1) , 0 )


# hydrogen mass fraction of the incoming gas
# 
X_H_MASSFRAC = 0.738

# preshock number density and velocity as function of the macrophysical parameters
# 
n_MdotMp(Mdot,Mp,Pop,ffill) = Mdot*1.893e20/(4*pi*( Rp_from_Mdot_Mp(Mdot,Mp,Pop)*7.15e9)**2.*ffill * -v_MdotMp(Mdot,Mp,Pop) * 1.673e-24/X_H_MASSFRAC)
v_MdotMp(Mdot,Mp,Pop) = StMoritz_v0simple(Mp, Rp_from_Mdot_Mp(Mdot,Mp,Pop))

# free-fall velocity
#   for free-fall from infinity ("simple": Racc = \infty)
# 
StMoritz_v0simple(Mp,Rp) = -1 * sqrt(2*6.67e-8*Mp*1.898e30/(Rp*7.15e9))


# ----------------------------------------------------------------------------------------

# simple test functions
# 
test_Rp_fit_Mdot_warm = "reset; set title 'WARM POPULATION'; set key top left Left rever title 'Mass (M_J)'; set aut; set log x; set xlabel 'dM/dt (M_E/yr)'; set ylabel 'Radius (R_J)'; plot [1e-5:1e-2] for [i in '20 15 10 5 3 1'] Rp_from_Mdot_Mp(x,i,'warm') w l lw 2 t i"

test_Rp_fit_Mp_warm = "reset; set title 'WARM POPULATION'; set key top left Left rever title 'dM/dt (M_E/yr)'; set aut; set xlabel 'Mass (M_J)'; set ylabel 'Radius (R_J)'; plot [0:30] for [i in '1e-1 2e-2 1e-2 1e-3 1e-4 1e-5'] Rp_from_Mdot_Mp(i,x,'warm') w l lw 2 t i"

test_Teff_fit_Mdot_warm = "reset; set title 'WARM POPULATION, f_{fill}=1'; ffill=1.0; set key top left Left rever title 'Mass (M_J)'; set aut; set log x; set xlabel 'dM/dt (M_E/yr)'; set ylabel 'T_{eff} (K)'; plot [1e-5:1e-2] for [i in '20 15 10 5 3 1'] Teff_from_Mdot_Mp(x,i,1000,'warm',ffill) w l lw 2 t i"

test_Teff_fit_Mp_warm = "reset; set title 'WARM POPULATION, f_{fill}=1'; ffill=1.0; set key top left Left rever title 'dM/dt (M_E/yr)'; set aut; set xlabel 'Mass (M_J)'; set ylabel 'T_{eff} (K)'; plot [0:30] for [i in '1e-1 2e-2 1e-2 1e-3 1e-4 1e-5'] Teff_from_Mdot_Mp(i,x,1000,'warm',ffill) w l lw 2 t i"

if(verbose>=1) print "# Loaded functions:"
if(verbose>=1) print "#   Rp_from_Mdot_Mp(Mdot,Mp,Pop)"
if(verbose>=1) print "#   Teff_from_Mdot_Mp(Mdot,Mp,Tint,Pop,ffill)"
if(verbose>=2) print "# and others (check for clash):"
if(verbose>=3) print "#   StMoritz_min()"
if(verbose>=3) print "#   StMoritz_max()"
if(verbose>=2) print "#   Tacc_MdotMp()"
if(verbose>=2) print "#   fdown1_Popsynthfit()"
if(verbose>=2) print "#   fdown_Popsynthfit()"
if(verbose>=2) print "#   n_MdotMp()"
if(verbose>=2) print "#   v_MdotMp()"
if(verbose>=3) print "#   StMoritz_v0simple()"
if(verbose>=2) print "# as well as the variables:"
if(verbose>=3) print "#   test_Rp_fit_Mdot_warm"
if(verbose>=3) print "#   test_Rp_fit_Mp_warm"
if(verbose>=3) print "#   test_Teff_fit_Mdot_warm"
if(verbose>=3) print "#   test_Teff_fit_Mp_warm"
if(verbose>=2) print "#   X_H_MASSFRAC"
if(verbose==2) print "# (a few other ones not shown)"
