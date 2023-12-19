;+
; NAME:
;   TARGET2BJD
; PURPOSE:
;   Converts a BJD in the target barycenter time to BJD in Barycentric
;   Dynamical Time (BJD_TDB).
;
; DESCRIPTION:
;   Corrects for the Roemer delay (light travel time) in
;   the target system (~30 seconds for typical Hot Jupiters). Most
;   importantly, this will naturally make primary, secondary, RV
;   observations, and phase variations self-consistent. 
;
;   Additionally, for a typical Hot Jupiter (3 hour transit in a 3 day
;   orbit), the ingress is delayed by ~0.15 seconds relative to the mid
;   transit time. For circular-orbits, the delay between mid transit
;   time and ingress is:
;
;   dt = a/c*(1-cos((sqrt((R_*+R_P)^2 - b^2)/a))
;
;   Falling for larger semi-major axis.
;
; INPUTS:
;   BJD_TARGET  - A scalar or array of BJDs in target time. Must be
;                 double precision.
;   INCLINATION - The inclination of the orbit
;   A           - The semi-major axis of the orbit (AU)
;   TP          - The time of periastron of the orbit (BJD_TARGET)
;   PERIOD      - The period of the orbit (days)
;   E           - Eccentricity of the orbit
;   OMEGA       - Argument of Periastron of the stellar orbit (radians)
;   
; OPTIONAL INPUTS:
;   Q          - The mass ratio of the targets (M1/M2). If not
;                specified, an infinite mass ratio is assumed (M1
;                stationary) (8 ms effect for Hot Jupiters).
;   C          - The speed of light, in AU/day. If not given,
;                initialized to 173.144483d0
;
; OPTIONAL KEYWORDS:
;   PRIMARY    - If set, the information comes from the position of
;                the primary (as in RV), and therefore the correction
;                will be the light travel time from the center of the
;                primary to the Barycenter -- analagous to the
;                difference between HJD and BJD in our solar system
;                (only ~8 ms for Hot Jupiters, but increasing with a).    
;                Otherwise, the correction will be the light
;                travel time from the smaller body to the barycenter
;                (as in transits) -- analagous to the difference
;                between JD and BJD in the Solar System.
;
;                NOTE: if Q is not specified and PRIMARY is, no
;                correction is applied.
;
; OUTPUTS:
;   BJD_TDB    - The time as it would flow in the Barycenter of the
;                solar system (BJD_TDB).
;
; LIMITATIONS:
;   We ignore the distance to the object (plane parallel waves), which
;   should have a similar effect as the distance plays in the BJD
;   correction (< 1 ms). We also ignore the systemic velocity, which
;   will compress/expand the period by a factor gamma/c.
;
; REVISION HISTORY:
; 2011/06: Written by Jason Eastman (OSU)

function target2bjd, bjd_target, inclination=inclination, a=a, tp=tp, $
                     period=period,e=e, omega=omega, q=q, primary=primary,c=c

if n_elements(q) eq 0 then q = !values.d_infinity

;; no correction necessary, already in the SSB frame
if  ~finite(q) and keyword_set(primary) then return, bjd_target

if n_elements(c) eq 0 then c = 173.144483d0 ;; AU/day

meananom = (2.d0*!dpi*(1.d0 + (bjd_target - Tp)/Period)) mod (2*!dpi) 

;; keplereq not vectorized if e has multiple elements
if n_elements(e) eq 1 then eccanom = exofast_keplereq(meananom, e) $
else begin
   nel = n_elements(meananom)
   if nel ne n_elements(e) then message, $
      "ERROR: e must be a scalar or be the same size as BJD_TARGET."
   eccanom = dblarr(nel)
   for i=0L, n_elements(meananom)-1 do $
      eccanom[i] = exofast_keplereq(meananom[i], e[i])
endelse
trueanom = 2d0*atan(sqrt((1d0 + e)/(1d0 - e))*tan(0.5d0*eccanom))

;; displacement from barycenter
if n_elements(q) ne 0 then begin
   if finite(q) then begin
      if keyword_set(primary) then factor = 1d0/(1d0+q) $ ;; a*factor = a1
      else factor = q/(1d0+q)                             ;; a*factor = a2
   endif else factor = 1d0
endif else factor = 1d0 ;; infinite mass ratio, a1=0, a2=a

;; distance from barycenter to target
r = a*(1d0-e^2)/(1d0+e*cos(trueanom))*factor

;; rotate orbit by omega
if ~keyword_set(primary) then om = omega + !dpi $
else om = omega

;; z grows with distance
z = r*sin(trueanom+om)*sin(inclination)

return, bjd_target - z/c

end    
    
