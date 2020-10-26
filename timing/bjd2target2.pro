;+
; NAME:
;   BJD2TARGET
; PURPOSE:
;   Converts a BJD in Barycentric Dynamical Time Time (BJD_TDB) to a BJD in
;   the target barycenter time.
;
; DESCRIPTION:
;   Iteratively corrects for the Roemer delay (light travel time) in
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
;   BJD_TDB     - A scalar or array of BJDs in TDB. Must be double
;                 precision.
;   INCLINATION - The inclination of the orbit
;   A           - The semi-major axis of the orbit (AU)
;   TP          - The time of periastron of the orbit (BJD_TARGET)
;   PERIOD      - The period of the orbit (days)
;   E           - Eccentricity of the orbit
;   OMEGA       - Argument of Periastron of the orbit (radians)
;   
; OPTIONAL INPUTS:
;   Q          - The mass ratio of the targets (M1/M2). If not
;                specified, an infinite mass ratio is assumed (M1
;                stationary) (8 ms effect for Hot Jupiters).
;   TOL        - The tolerance, in days. The iterative procedure will
;                stop after the worst-case agreement is better than this.
;                Default = 1d-8 (1 ms).
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
;   BJD_TARGET - The time as it would flow in the Barycenter of the target.
;
; LIMITATIONS:
;   We ignore the distance to the object (plane parallel waves), which
;   should have a similar effect as the distance plays in the BJD
;   correction (< 1 ms). We also ignore the systemic velocity, which
;   will compress/expand the period by a factor gamma/c.
;
; REVISION HISTORY:
; 2011/06: Written by Jason Eastman (OSU)

function bjd2target2, bjd_tdb, inclination=inclination, a=a, tp=tp, $
                      period=period,e=e, omega=omega, q=q, tol=tol, $
                      primary=primary

;; short circuit this (expensive) calculation.
;return, bjd_tdb

if n_elements(tol) eq 0 then tol = 1d-8 ;; 1 ms

c = 173.144483d0 ;; AU/day
bjd_target = bjd_tdb

depthold = 0

;; iterative process to find the BJD_Target that corresponds to the BJD_TDB
;; completes in ~3 iterations
repeat begin

   target_new = target2bjd(bjd_target,inclination=inclination, a=a, tp=tp, $
                           period=period,e=e,omega=omega,q=q,primary=primary)

   diff = bjd_tdb-target_new
   bjd_target += diff

endrep until max(abs(diff)) lt tol

return, bjd_target

end    
    
