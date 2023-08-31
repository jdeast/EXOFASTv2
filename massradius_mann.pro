;+
; NAME:
;   MASSRADIUS_MANN
;
; PURPOSE: 
;   Uses the Mann relation
;   (https://ui.adsabs.harvard.edu/abs/2015ApJ...804...64M/abstract)
;   the stellar mass and radius given Ks and [Fe/H].
;
; CALLING SEQUENCE:
;    massradius_mann, ks0, mstar, rstar, feh=feh,
;                     sigma_rstar=sigma_rstar, sigma_mstar=sigma_mstar, 
;                     distance=distance,mann19=mann19
;
; INPUTS:
;    KS0  - The K-band magnitude. If distance is given, this is
;           interpreted as the apparent magnitude. If distance is not
;           given, it's interpreted as the absolute magnitude.
; OPTIONAL INPUTS
;    FEH  - The stellar metallicity. If supplied, we use equation
;           5. If not, we use equation 4.
;    DISTANCE - The stellar distance (pc). If supplied, we
;               assume KS0 is the apparent magnitude and we compute the
;               absolute magnitude
; OPTIONAL KEYWORDS - 
;    MANN19 - If set, uses the Mann 2019 relations to estimate the mass.
; OUTPUTS:
;    MSTAR       - The stellar mass, in solar masses
;    RSTAR       - The stellar radius, in solar radii
;    SIGMA_MSTAR - The stellar mass uncertainty, in solar masses
;    SIGMA_RSTAR - The stellar radius uncertainty, in solar radii
;
; MODIFICATION HISTORY
; 
;  2023/08 -- Public release -- Jason Eastman
;-
pro massradius_mann, ks0, mstar, rstar, feh=feh, $
                     sigma_rstar=sigma_rstar, sigma_mstar=sigma_mstar,$
                     distance=distance,mann19=mann19

;; compute the absolute magnitude if given the apparent magnitude
if n_elements(distance) eq 0 then begin
   ks = ks0
endif else begin
   ks = ks0 - 2.5d0*alog10((distance/10d0)^2)
endelse

if n_elements(feh) eq 0 then begin
   ;; eq 4, table 1
   ai = [1.9515d0, -0.3520d0, 0.01680d0]            
   rstar = ai[0] + ai[1]*ks + ai[2]*ks^2
   sigma_rstar = rstar*0.0289d0
endif else begin
   ;; eq 5, table 1
   ai = [1.9305d0, -0.3466d0, 0.01647d0,0.04458d0] 
   rstar = ai[0] + ai[1]*ks + ai[2]*ks^2 + ai[3]*feh
   sigma_rstar = rstar*0.027d0
endelse

if keyword_set(mann19) then begin
   zp = 7.5d0
   if n_elements(feh) eq 0 then begin
      ;; eq 4, table 6 (n=5)
      bi = [-0.642d0,-0.208d0,-8.43d-4,7.87d-3,1.42d-4,-2.13d-4]
      mstar = 10d0^(bi[0] + bi[1]*(ks-zp) + bi[2]*(ks-zp)^2 + bi[3]*(ks-zp)^3 + bi[4]*(ks-zp)^4 + bi[4]*(ks-zp)^5)
      sigma_mstar = mstar*0.020d0
   endif else begin
      ;; eq 5, table 6 (n=5)
      bi = [-0.647d0,-0.207d0,-6.53d-4,7.13d-3,1.84d-4,-1.6d-4,-0.0035]
      mstar = (1d0+feh*bi[6])*10d0^(bi[0] + bi[1]*(ks-zp) + bi[2]*(ks-zp)^2 + bi[3]*(ks-zp)^3 + bi[4]*(ks-zp)^4 + bi[4]*(ks-zp)^5)
      sigma_mstar = mstar*0.021d0
   endelse
endif else begin
   ;; Mann 2015 et al., eq 10, table 1
   bi = [0.5858d0, 0.3872d0, -0.1217d0, 0.0106d0, -2.7262d-4] 
   mstar = bi[0] + bi[1]*ks + bi[2]*ks^2 + bi[3]*ks^3 + bi[4]*ks^4
   sigma_mstar = mstar*0.018d0
endelse



end
