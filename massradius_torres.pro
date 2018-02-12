;+
; NAME:
;   MASSRADIUS_TORRES
;
; PURPOSE: 
;   Uses the Torres relation
;   (http://adsabs.harvard.edu/abs/2010A%26ARv..18...67T) to determine
;   the stellar mass and radius given the logg, Teff, and [Fe/H].
;
;   NOTE the errors in the empirical model of ulogm = 0.027d0, ulogr = 0.014d0
;
; CALLING SEQUENCE:
;    massradius_torres, logg, teff, feh, mstar, rstar
;
; INPUTS:
;    LOGG - The log of the stellar surface gravity
;    TEFF - The stellar effective temperature
;    FEH  - The stellar metalicity
;
; OUTPUTS:
;    MSTAR - The stellar mass, in solar masses
;    RSTAR - The stellar radius, in solar radii
;
; MODIFICATION HISTORY
; 
;  2012/06 -- Public release -- Jason Eastman (LCOGT)
;-
pro massradius_torres, logg, teff, feh, mstar, rstar

;; coefficients from Torres, 2010
ai = [1.5689d,1.3787d,0.4243d,1.139d,-0.14250d,0.01969d,0.10100d]
bi = [2.4427d,0.6679d,0.1771d,0.705d,-0.21415d,0.02306d,0.04173d]
;ulogm = 0.027d0
;ulogr = 0.014d0

X = alog10(teff) - 4.1d0
   
sz = size(teff)
if sz[0] eq 1 then begin
   coeffs = transpose([[dblarr(sz[1])+1.d0],[X],[X^2],[X^3],[logg^2],$
                       [logg^3d0],[feh]])
   logm = total((replicate(1d0,sz[1])##ai)*coeffs,1)
   logr = total((replicate(1d0,sz[1])##bi)*coeffs,1)
   nran = sz[1]
endif else if sz[0] eq 0 then begin
   coeffs = [1.d0,X,X^2d0,X^3,logg^2d0,logg^3d0,feh]
   logm = total(ai*coeffs)
   logr = total(bi*coeffs)
   nran = 1
endif else message, 'Teff, Fe/H, logg  must be scalars or 1D vectors'

mstar = 10d0^logm
rstar = 10d0^logr

end
