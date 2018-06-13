;+
; NAME:
;   BJD2TARGET
; PURPOSE:
;   Iteratively calls TARGET2BJD to convert a BJD in Barycentric
;   Dynamical Time (BJD_TDB) to a BJD in the target barycenter
;   time (BJD_TARGET) within TOL days.
;
; DESCRIPTION:
;   The opposite of TARGET2BJD; see description there.
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
;                specified, an infinite mass ratio is assumed (M1 is
;                stationary at the barycenter) (8 ms effect for Hot
;                Jupiters).
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
;-

function bjd2target, bjd_tdb, inclination=inclination, a=a, tp=tp, $
                     period=period,e=e, omega=omega, q=q, tol=tol, $
                     primary=primary,pars=pars, c=c

if n_elements(tol) eq 0 then tol = 1d-8 ;; 1 ms
bjd_target = bjd_tdb

niter = 0
;; iterative process to find the BJD_Target that corresponds to the BJD_TDB
;; completes in ~3 iterations
repeat begin

   target_new = target2bjd(bjd_target,inclination=inclination, a=a, tp=tp, $
                           period=period,e=e,omega=omega,q=q,primary=primary,c=c)

   diff = bjd_tdb-target_new
   bjd_target += diff
   niter++
   
;   if ~finite(max(abs(diff))) then stop

   if niter gt 100 then message, 'Not converging; this is a rare bug usually associated with poorly constrained parameters. Try again or consider imposing some priors for poorly constrained parameters'

endrep until max(abs(diff)) lt tol

return, bjd_target

end    
    
