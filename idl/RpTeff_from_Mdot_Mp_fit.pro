;  AMDG
; JMJ-V!
; FdS 24.01.2020
; 
; This file contains two functions giving the radius and the Teff of forming
; planets. They are defined in Aoyama et al. (2020)
; 
; A.)
; The function Rp_from_Mdot_Mp(Mdot,Mp,Pop) returns the fit of the radius as a 
; function of accretion rate Mdot and planet mass Mp for the cold- (Pop='cold')
; or warm-start (Pop='warm') populations of Mordasini et al. (2012). Plotted in
; Fig. 1 of Aoyama et al. (2020).
; 
; B.)
; The function Teff_from_Mdot_Mp(Mdot,Mp,Tint,Pop,ffill) returns the fit of
; of the Teff of the accreting region of the planet's surface. It ensures
; approximate energy conservation. The details of the fit might change in 
; future iterations of the model. See §2.2.2 and Figs. 2 and 3.
; 
; The non-accreting region is in principle emitting at Teff = Tint
; but this is not relevant for this code.
; 
; The fits are valid for approximately (please see the paper):
; 
;          1 MJ <  Mp  < 30 MJ
;    1e-5 ME/yr < Mdot < 2e-2 ME/yr
; 
; -- Input --
;   Mdot: in ME/yr = 1.893e20 g/s
;     Mp: in MJ = 1.898e30 g
;    Pop: 'warm' or 'cold'
; 
; -- Output --
;     Rp: in RJ = 7.15e9 cm
;   Teff: in K
; 
; To test the functions, you can do the following:
; 
; IDL> .com Rp_from_Mdot_Mp_Popsynthfit.pro
; IDL> test_Rp_fit, mode='Mdot', Pop='warm'
; 
; and so on, varying mode to 'Mp' and Pop to 'cold',
; or simply
; 
; IDL> .run Rp_from_Mdot_Mp_Popsynthfit.pro
; 
; If you use these functions, please cite:
; 
;   Aoyama, Marleau, Mordasini & Ikoma (subm.)
;   Spectral appearance of the planetary-surface accretion shock:
;     Global spectra and hydrogen-line profiles and luminosities
;   https://arxiv.org/abs/2011.06608
; 
; Comments and requests welcome!
; 
;   Gabriel-Dominique Marleau (c) 2021
;   @uni-tuebingen.de with gabriel.marleau in front
; 
; -------------------------------------------------------------
; 
;  Note: use real() to be able to use the function at Mp < 1
;  for esthetic reasons
; 
function Rp_from_Mdot_Mp, Mdot, Mp, Pop
    
    lm2 = alog10(Mdot/1e-2)
    if (Pop eq 'warm') then begin
        Rp = real_part( 0.411-0.244*lm2+3.45*exp(0.762*lm2) + (-0.489-0.0961*lm2+0.652*exp(0.353*lm2))*(Mp-1.0) + (-0.228-0.00106*lm2+0.226*exp(0.000220*lm2))*(Mp-1.0)^2. )
    endif else if (Pop eq 'cold') then begin
        Rp = real_part( 1.53+0.111*lm2+1.06*exp(0.906*lm2) + (-0.195-0.0307*lm2+0.0977*exp(0.000695*lm2))*(Mp-1.0) + (-0.250+0.000276*lm2+0.254*exp(0.000214*lm2))*(Mp-1.0)^2. )
    endif else begin
        print,' * Pop should be one of "warm" or "cold" *'
        stop
    endelse
    
    return, Rp
     
end

function Teff_from_Mdot_Mp, Mdot, Mp, Tint, Pop, ffill

    return, (Tint^4. +              fdown_Popsynthfit(alog10(n_MdotMp(Mdot,Mp,Pop,ffill)),v_MdotMp(Mdot,Mp,Pop)/(-1e5)) * Tacc_MdotMp(Mdot,Mp,Pop,ffill)^4.)^0.25

end

; -----------------------------------------------------------------------------

; 
; the "accretion temperature" of the accreting region
;   - Here it is correct to divide by sigma_SB=5.67e-5 cgs and not 4*sigma_SB = a*c
;     because we want the effective temperature, not the gas temperature
;     (see Marleau et al. 2019; http://adsabs.harvard.edu/abs/2019ApJ...881..144M)
;   - The rest of the planet's surface is emitting only with Tint
; 
function Tacc_MdotMp, Mdot,Mp,Pop,ffill
    return, ( 6.67d-8*Mp*1.898d30*Mdot*1.893d20/ (ffill * 4*!dpi*(Rp_from_Mdot_Mp(Mdot,Mp,Pop)*7.15d9)^3. * 5.67d-5) )^0.25d0
end

; 
; fit to the data of Aoyama et al. (2018) and Aoyama et al. (2020; https://arxiv.org/abs/2011.06608)
; 
function fdown1_Popsynthfit, lgn0, v0

    return, 0.703752d0-0.0967987d0*(lgn0-12)-0.0254579d0*abs(lgn0-12)^2.  +  (-0.00527886d0-0.00146833d0*(lgn0-12)-0.000321504d0*abs(lgn0-12)^2.)*(v0-100)  -9.91492d-06*abs(v0-100)^2.

end

; limit fdown to [0,1]
; 
function fdown_Popsynthfit, lgn0, v0
    return, ( fdown1_Popsynthfit(lgn0,v0) < 1 ) > 0
end

; preshock number density and velocity as function of the macrophysical parameters
; 
function n_MdotMp, Mdot, Mp, Pop,ffill
    
    ; hydrogen mass fraction of the incoming gas
    ; 
    X_H_MASSFRAC = 0.738
    
    return, Mdot*1.893d20/(4*!dpi*( Rp_from_Mdot_Mp(Mdot,Mp,Pop)*7.15d9)^2.*ffill * (-v_MdotMp(Mdot,Mp,Pop)) * 1.673e-24/X_H_MASSFRAC)
end

function v_MdotMp, Mdot, Mp, Pop
    return, StMoritz_v0simple(Mp, Rp_from_Mdot_Mp(Mdot,Mp,Pop))
end

; free-fall velocity
;   for free-fall from infinity ("simple": Racc = \infty)
; 
function StMoritz_v0simple, Mp, Rp
    return, -1 * sqrt(2*6.67d-8*Mp*1.898d30/(Rp*7.15d9))
end

; -----------------------------------------------------------------------------


; test routine for the radius fit
; 
pro test_Rp_fit, mode, Pop
    th=1.8
    cs=2
    
    ; plot with Mdot on the x axis
    if (mode eq 'Mdot') then begin
        MM = [20, 15, 10, 5, 3, 1]
        Mdot = 10^(-5. + (-1+5)/50.*dindgen(51))
        
        ; set up plot
        plot, [1],[1], /nodata, /xst,/yst, /xlog, xrange=[1e-5,1e-1], /ylog, yrange=[1,8], yticks=3,ytickv=[1,2,4,8], $
            xtit=textoidl('!6dM/dt (M_E/yr)'), ytit=textoidl('Radius (R_J)'),xth=th,yth=th, chars=cs, tit=Pop
        
        print, ' Masses (M_J):'
        print, MM
        
        ; plot each mass
        for i=0,n_elements(MM)-1 do begin
            Rp = Rp_from_Mdot_Mp(Mdot, MM[i], Pop)
            oplot, Mdot, Rp, thick=th
        endfor
    endif
    
    ; plot with Mp on the x axis
    if (mode eq 'Mp') then begin
        Mdot = [2e-2, 1e-2, 1e-3, 1e-4, 1e-5]
        MM = 0.5 + (30.-0.5)/50.*dindgen(51)

        ; set up plot
        plot, [1],[1], /nodata, /xst,/yst, xrange=[0,30], /ylog, yrange=[1,8], yticks=3,ytickv=[1,2,4,8], $
            xtit=textoidl('!6Mass (M_J)'), ytit=textoidl('Radius (R_J)'),xth=th,yth=th, chars=cs, tit=Pop
          
        print, ' Accretion rates (M_E/yr):'
        print, Mdot
        
        ; plot each accretion rate
        for i=0,n_elements(Mdot)-1 do begin
            Rp = Rp_from_Mdot_Mp(Mdot[i], MM, Pop)
            oplot, MM, Rp, thick=th
        endfor
    endif
end

; -----------------------------------------------------------------------------

; test routine for the Teff fit
; 
pro test_Teff_fit, mode, Pop, ffill
    th=1.8
    cs=2
    Title = textoidl('T_{eff} of hot spot. Population: ')
    
    ; effectively a minimum temperature
    Tint = 1000
    
    ; plot with Mdot on the x axis
    if (mode eq 'Mdot') then begin
        MM = [20, 15, 10, 5, 3, 1]
        Mdot = 10^(-5. + (-1+5)/50.*dindgen(51))
        
        ; set up plot
        plot, [1],[1], /nodata, /xst,/yst, /xlog, xrange=[1e-5,1e-1], /ylog, $
        yrange=[700,7000], $ ; , yticks=3 , ytickv=[1,2,4,8], $
            xtit=textoidl('!6dM/dt (M_E/yr)'), ytit=textoidl('T_{eff} (K)'),xth=th,yth=th, chars=cs, tit=Title+Pop
        
        print, ' Masses (M_J):'
        print, MM
        
        ; plot each mass
        for i=0,n_elements(MM)-1 do begin
            Teff = Teff_from_Mdot_Mp(Mdot, MM[i], Tint, Pop, ffill)
            oplot, Mdot, Teff, thick=th
        endfor
        oplot, Mdot, Mdot*0 + Tint, lines=1, col=254, thick=th
    endif
    
    ; plot with Mp on the x axis
    if (mode eq 'Mp') then begin
        Mdot = [2e-2, 1e-2, 1e-3, 1e-4, 1e-5]
        MM = 0.5 + (30.-0.5)/50.*dindgen(51)
        
        print, ' Accretion rates (M_E/yr):'
        print, Mdot
        
        ; set up plot
        plot, [1],[1], /nodata, /xst,/yst, xrange=[0,30], /ylog, $
        yrange=[500,1e4], $
            xtit=textoidl('!6Mass (M_J)'), ytit=textoidl('T_{eff} (K)'),xth=th,yth=th, chars=cs, tit=Pop
          
        ; plot each accretion rate
        for i=0,n_elements(Mdot)-1 do begin
            Teff = Teff_from_Mdot_Mp(Mdot[i], MM, Tint, Pop, ffill)
            oplot, MM, Teff, thick=th
        endfor
        oplot, MM, MM*0 + Tint, lines=1, col=254, thick=th
    endif
end

; -----------------------------------------------------------------------------

; test routines for internal functions

pro test_fdown
    th=1.8
    
    lgn0 = 14
    v0 = 20 + 180*dindgen(101)/100.
    
    fd1 = fdown1_Popsynthfit(lgn0, v0)
    fd  = fdown_Popsynthfit(lgn0, v0)
    
    plot, v0, fd1, $
        xtit=textoidl('!6v_0 (km/s)'), ytit=textoidl('f_{down}'), $
        lines=0, $
        xth=th,yth=th, chars=2, thick=th
    oplot,v0, fd, lines=1, col=254, thick=th
end


; =============================================================================

; to run the tests:
; IDL> .run Rp_from_Mdot_Mp_Popsynthfit.pro
; IDL> .cont
; IDL> .cont

test_fdown
stop

Modus = 'Mdot'  ; which variable is the independent one
Pop = 'cold'
ffill = 1.0

test_Rp_fit, Modus, Pop
stop
test_Teff_fit, Modus, Pop,ffill

end
