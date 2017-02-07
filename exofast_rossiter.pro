pro exofast_rossiter, x, y, u1, p, vsini, lambda, deltarv, z=z
;+
; NAME:
;   EXOFAST_ROSSITER
;
; PURPOSE: 
;   This routine implements the Ohta et al 2005 approximation for the
;   Rossiter-Mclaughlin effect, which assumes linear limb darkening
;   and is accurate to ~1 m/s for slowly rotating stars (see Hirano 2011).
;
; CALLING SEQUENCE:
;   exofast_rossiter, x, y, u1, p, vsini, lambda, deltarv, z=z
;
; INPUTS:
;    x      - array of star-planet distances in the x direction
;             (stellar radii).
;    y      - array or scalar of star-planet distance(s) in the y
;             direction. ie, the impact parameter (stellar radii)
;    u1     - linear limb-darkening coefficient -- epsilon in paper
;    p      - occulting star size (stellar radii) -- gamma in paper
;    vsini  - projected stellar surface velocity (same as output units)
;    lambda - projected spin-orbit angle (radians)
;    
; OPTIONAL INPUTS:
;    z      - array of star-planet distances in the z direction. 
;               -- Positive => in front of star (RM effect)
;               -- Negative => behind star (no effect)
;               -- assumed positive if not specified
; OUTPUTS:
;    deltarv - Approximate Change in RV at each z due to the Rossiter
;              Mclaughlin effect (units same as vsini) 
;                     *** accurate to ~1 m/s ***
;
; MODIFICATION HISTORY
;  2009/04/30 -- Jason Eastman (Ohio State University)
;     Loose translation from Scott Gaudi's (OSU) fortran code
;  2012/08 -- Jason Easmtan (LCOGT)
;     Make sign convention ahere to standard (now consistent with EXOFAST_GETB).
;-

nx = n_elements(x)
if n_elements(z) eq 0 then z = dblarr(nx)+1d0
deltarv = dblarr(nx)
rho = sqrt(x^2+y^2)

xp = x*cos(lambda)-y*sin(lambda)
;yp = x*sin(lambda)+y*cos(lambda)
;theta = dindgen(1000)/1000d0*2d0*!dpi
;plot, cos(theta),sin(theta),/iso, xrange=[-2,2],yrange=[-2,2]
;oplot, xp, yp

etap = rho-1.d0         ;; equation 28

;; transit phase
transit = where(rho le (1d0-p) and z gt 0d0)
if transit[0] ne -1 then begin
    w1 = 0d0                    ;; equation 43
    w2 = sqrt(1-rho[transit]^2) ;; equation 44
    ;; equation 40
    deltarv[transit]=-vsini*xp[transit]*p^2*(1d0-u1*(1d0-w2))/$
      (1d0-p^2-u1*(1d0/3d0-p^2*(1d0-w1)))
endif

;; ingress/egress phase
inegress = where(abs(etap) lt p and z gt 0)
if inegress[0] ne -1 then begin
    ;; equation 33
    x0 = 1.d0 - (p^2-etap[inegress]^2)/(2d0*(1d0+etap[inegress])) ;; equation 33
    zeta = 1.d0 + etap[inegress] - x0 ;; equation 34
    z0 = sqrt(1-x0^2)       ;; equation 33
    xc=x0+(zeta-p)/2d0      ;; equation 50
    w2 = sqrt(2d0*p-p^2)    ;; equation 44
    w3 = 0.d0               ;; equation 48

    ;; equation 51
    g1 = (1d0-xc^2)*asin(sqrt((p^2-(xc-1-etap[inegress])^2)/(1.d0-xc^2))) + $
      sqrt((p^2-(xc-1.d0-etap[inegress])^2)*(1.d0-xc^2-p^2+$
                                             (xc-1.d0-etap[inegress])^2))
    g2 = sqrt(2.d0*p^3*(1.d0-p)) + p*(2.d0-p)*asin(sqrt(p/(2.d0-p)))

    w4 = (!dpi/2.d0)*p*(p-zeta)*xc*g1/g2*w2 ;; equation 49

    ;; equation 45
    num=(1.d0-u1)*(-z0*zeta+p^2*acos(zeta/p))+$
      u1*w4/(1.d0+etap[inegress])

    den1 = !dpi*(1.d0-1.d0/3.d0*u1)-u1*w3
    den2 = -(1.d0-u1)*(asin(z0)-$
               (1.d0+etap[inegress])*z0+p^2*acos(zeta/p))
    deltarv[inegress]=-vsini*xp[inegress]*num/(den1+den2)
endif

end

