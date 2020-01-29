;  AMDG
; JMJ-V!
; FdS 24.01.2020
; 
; The function Rp_from_Mdot_Mp_Popsynthfit(Mdot,Mp,warmPop) returns the fit,
; given in Aoyama, Marleau, Mordasini & Ikoma (2020),
; to the radius as a function of accretion rate Mdot and planet mass Mp
; for the cold- (Pop='cold') or warm-start (Pop='warm') population
; of Mordasini et al. (2012).
; The fit is valid for approximately (please see the paper):
; 
;          1 MJ <  Mp  < 30 MJ
;    1e-5 ME/yr < Mdot < 2e-2 ME/yr
; 
; -- Input --
;   Mdot: in ME/yr = 1.893e20 g/s
;     Mp: in MJ = 1.898e30 g
; 
; -- Output --
;     RP: in RJ = 7.15e9 cm
; 
; To test the function, you can do the following:
; 
; .com Rp_from_Mdot_Mp_Popsynthfit.gnu
; 
; test_Rp_fit, mode='Mdot',Pop='warm'
; 
; If you use this function, please cite: Aoyama et al. (2020). Thank you!
; 
; Comments and requests welcome!
; 
; Gabriel-Dominique Marleau (c) February 2020
; gabriel.marleau@uni-tuebingen.de, gabriel.marleau@space.unibe.ch
; 
; -------------------------------------------------------------
; 
;  Note: use real() to be able to use the function at Mp < 1
;  for esthetic reasons
; 
function Rp_from_Mdot_Mp_Popsynthfit, Mdot, Mp, Pop
    
    lm2 = alog10(Mdot/1e-2)
    if (Pop eq 'warm') then begin
        Rp = real_part( 0.411-0.244*lm2+3.45*exp(0.762*lm2) + (-0.489-0.0961*lm2+0.652*exp(0.353*lm2))*(Mp-1.0) + (-0.228-0.00106*lm2+0.226*exp(0.000220*lm2))*(Mp-1.0)^2. )
    endif else begin
        Rp = real_part( 1.53+0.111*lm2+1.06*exp(0.906*lm2) + (-0.195-0.0307*lm2+0.0977*exp(0.000695*lm2))*(Mp-1.0) + (-0.250+0.000276*lm2+0.254*exp(0.000214*lm2))*(Mp-1.0)^2. )
    endelse
    
    return, Rp
     
end

; simple test routine
; 
pro test_Rp_fit, mode, Pop
    th=1.7
    cs=2
    
    ; plot with Mdot on the x axis
    if (mode eq 'Mdot') then begin
        MM = [20, 15, 10, 5, 3, 1]
        Mdot = 10^(-5. + (-1+5)/50.*dindgen(51))
        
        ; set up plot
        plot, [1],[1], /nodata, /xst,/yst, /xlog, xrange=[1e-5,1e-1], /ylog, yrange=[1,8], yticks=3,ytickv=[1,2,4,8], $
            xtit=textoidl('!6dM/dt (M_E/yr)'), ytit=textoidl('Radius (R_J)'),xth=th,yth=th, chars=cs, tit=Pop
        
        ; plot each mass
        for i=0,n_elements(MM)-1 do begin
            Rp = Rp_from_Mdot_Mp_Popsynthfit(Mdot, MM[i], Pop)
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
          
        ; plot each accretion rate
        for i=0,n_elements(Mdot)-1 do begin
            Rp = Rp_from_Mdot_Mp_Popsynthfit(Mdot[i], MM, Pop)
            oplot, MM, Rp, thick=th
        endfor
    endif
end
