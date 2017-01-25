function exofast_rv, bjd, TPeriastron, period, V0, K, e0, omega0, $
                     slope=slope,t0=t0,rossiter=rossiter,i=i,a=a,u1=u1,p=p, $
                     vsini=vsini, lambda=lambda,deltarv=deltarv
;+
; NAME:
;   exofast_rv
;
; PURPOSE: 
;   This function returns the radial velocity at each BJD. Optionally
;   includes the Rossiter McLaughlin effect, using the Ohta et al 2005
;   approximation good to ~1 m/s.
; 
; CALLING SEQUENCE:
;   result = exofast_rv(bjd, TPeriastron, period, V0, K, e, omega, $
;                       i=i, a=a, u1=u1, p=p, vsini=vsini, lambda=lambda, $
;                       r_star=r_star, slope=slope)
;
; INPUTS:
;   bjd         - Barycentric Julians dates for desired impact parameters
;   tperiastron - periastron passage time (BJD)
;   Period      - Period of orbit (days)
;   V0          - systemic velocity (m/s)
;   K           - radial velocity semi-amplitude (m/s)
;
; OPTIONAL INPUTS:  
;   e           - eccentricity of orbit (0 if not specified)
;   omega       - orbit's argument of periastron (radians) 
;                 - required if e is specified
;                 - assumed to be pi/2 if e not specified
;   slope       - the slope of the RV signal
;   t0          - the time at which V0 is referenced, if a slope is
;                 given. If not given, (maxtime + mintime)/2 is
;                 assumed.
;                 NOTE: if comparing two data sets with a slope and
;                 different times, this must be specified.
;
; OPTIONAL KEYWORDS:
;   rossiter    - If set, calculates the RM effect. All values below
;                 must also be specified.
;
; Optional parameters (required for Rossiter-McLaughlin effect)
;   i           - inclination of orbit (radians)
;   a           - semi-major axis (R_star)
;   u1          - linear limb-darkening coefficient
;   p           - occulting star size (stellar radii)
;   vsini       - projected stellar surface velocity (m/s)
;   lambda      - projected spin-orbit angle (radians)   
;
; OUTPUTS:
;    result     - the radial velocity of the star at each point
;
; MODIFICATION HISTORY 
;  2009/05/01 -- Jason Eastman (Ohio State University)
;  2013/10/15 -- Clarified documentation
;-

;; calculate the mean anomaly corresponding to each observed time
meananom = 2.d0*!dpi*(1.d0 + (bjd - Tperiastron)/Period mod 1)

;; if eccentricity is given, integrate the orbit
if n_elements(e0) ne 0 then begin
    if n_elements(omega0) eq 0 then message, $
      'ERROR: omega must be specified if e is specified'

    if e0 lt 0 then begin
        e = -e0
        omega = omega0 + !dpi
    endif else begin
        e = e0
        omega = omega0
    endelse

    eccanom = exofast_keplereq(meananom, e)
    trueanom = 2.d0*atan(sqrt((1.d0 + e)/(1.d0 - e))*tan(eccanom/2.d0))
endif else begin
    trueanom = meananom 
    e=0.d0
endelse

;; standard definition of omega for circular orbits
if n_elements(omega) eq 0 then omega = !dpi/2.d0

;; RV signal with no RM effect
rv = K*(cos(trueanom+omega) + e*cos(omega)) + V0

;; add a slope if desired
if n_elements(slope) ne 0 then begin
    if n_elements(t0) eq 0 then begin
        mintime = min(bjd,max=maxtime)
        t0 = (maxtime+mintime)/2.d0
    endif
    rv += (bjd - t0)*slope
endif

;; Calculate the RM effect
if keyword_set(rossiter) then begin
    if n_elements(i) eq 0 or n_elements(a) eq 0 or n_elements(u1) eq 0 or $
      n_elements(p) eq 0 or n_elements(vsini) eq 0 or n_elements(lambda) eq 0 $
      then message, 'ERROR: a, i, u1, p, vsini, and lambda must be ' + $
      'specified in order to calculate the Rossiter McLaughlin effect'

    ;; calculate the corresponding (x,y) coordinates of planet
    r = a*(1d0-e^2)/(1d0+e*cos(trueanom))
    
    ;; as seen from observer
    x = -r*cos(trueanom + omega)
    tmp = r*sin(trueanom + omega)
    y =  -tmp*cos(i)
    z =  tmp*sin(i)

    exofast_rossiter, x, y, u1, p, vsini, lambda, deltarv, z=z
    rv += deltarv

endif

return, rv

end
