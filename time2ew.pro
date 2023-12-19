; NAME:
;   TIME2EW
;
; PURPOSE:
;   Compute e and omega from the transit times and durations using
;   Kopal 1946 equations.
;
; INPUTS:
;   TC - The time of primary transit (BJD_TDB)
;
;   TS - The time of secondary eclipse (BJD_TDB)
;
;   PERIOD - The period of the orbit (days)
;
; OPTIONAL INPUTS
;   
;   ESINWIN - e*sin(omega). If not given, computed from
;             transit/eclipse durations.
;   
;   TFWHMP  - The primary transit duration from 1.5 contact to 3.5
;             contact (days). If not given, esinw=0.
;
;   TFWHMS  - The secondary eclipse duration from 1.5 contact to 3.5
;             contact (days). If not given, esinw=0.
;
;   P       - Rp/Rstar. Used for a second order correction in
;             esinw. If not given, no second order correction to esinw.
;
;   INC     - Inclination (radians). Used for a second order correction in
;             esinw. If not given, no second order correction to esinw.
;   
;   TOL     - Maximum error in the numerical solution for
;             ecosw. Default is 1d-9
;  
;   MAXITER - Maximum number of iterations before it throws an
;             error. Default is 100.
;
; OUTPUTS:
;   
;   E       - The orbital eccentricity
;  
;   OMEGA   - Argument of periastron of the star's orbit, in radians
; 
;   ECOSW   - e*cos(omega_*)
;
;   ESINWOUT- e*sin(omega_*)
;
; MODIFICATION HISTORY:
;   2023/12/18 - Jason Eastman (CfA), Written
;-
function time2ew, tc, ts, period, esinwin=esinwin, $
                  tfwhmp=tfwhmp, tfwhms=tfwhms, p=p, inc=inc, $
                  e=e, omega=omega, ecosw=ecosw, esinwout=esinwout,$
                  tol=tol, maxiter=maxiter

;; Kopal 1946 eq 133
if n_elements(esinwin) eq 0 then begin
   ;; determine (approximate) esinw from transit/eclipse durations
   if n_elements(tfwhmp) ne 0 and n_elements(tfwhms) ne 0 then begin  
      if n_elements(p) ne 0 and n_elements(inc) ne 0 then begin
         num = 1d0 - (1d0/(1d0+p))^2*(cos(inc))^2
         den = 1d0 -   (p/(1d0+p))^2*(cos(inc))^2
         beta = sqrt(num/den)
      endif else beta = 1d0
      esinw = (tfwhms*beta - tfwhmp)/(tfwhms*beta + tfwhmp)
   endif else esinw = 0d0
endif else esinw = esinwin

;; Kopal says this could be done as a parabolic case of Kepler's
;; equation, but my solver can't handle parabolic orbits (e=1)
;; and this doesn't recover full range of omega
;psi = exofast_keplereq(2d0*!dpi/period*(ts-tc),1d0)
;tanx2 = tan((psi-!dpi)/2d0)^2 
;e = sqrt((tanx2 + esinw^2)/(1d0+tanx2))
;omega = arcsin(esinw/e)

;; numerically solve for ecosw
minecosw = -1d0
maxecosw = 1d0

if n_elements(tol) eq 0 then tol=1d-9
if n_elements(maxiter) eq 0 then maxiter = 100
iter = 0
repeat begin

   if iter gt maxiter then begin
      print, 'not converging on value for ecosw'
      stop
   endif

   ecosw = (minecosw+maxecosw)/2d0
   e = sqrt(ecosw^2+esinw^2)
   omega = atan(esinw,ecosw)

   ;; Kopal, 1946 eq 114
   psi = !dpi + 2d0*atan(ecosw/(sqrt(1d0-e^2)))

   ;; Kopal 1946, eq 115
   residual = 2d0*!dpi/period*(ts-tc) - psi + sin(psi)

   if residual gt 0 then minecosw = ecosw $
   else maxecosw = ecosw

   iter++

endrep until abs(residual) lt tol
esinwout=esinw

return, ecosw

end
