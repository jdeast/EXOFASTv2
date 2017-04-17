;+
; NAME:
;   EXOFAST_GETB
;
; PURPOSE: 
;   This function returns the impact parameter as a function of BJD
;   given orbital elements of the planet, assuming a keplerian orbit.
;   Optionally, the 3 space coordinates can be returned.
; 
; CALLING SEQUENCE:
;   result = getb(jds, i=i,a=a,tperiastron=tperiastron,period=p,$
;                 e=e,omega=omega, x=x, y=y, z=z)
;
; INPUTS:
;   bjd         - Barycentric Julians dates for desired impact
;                 parameters (ideally in the target's
;                 barycentric frame; see BJD2TARGET)
;   i           - inclination of orbit (radians)
;   a           - semi-major axis, a = a1+a2 (in units of R_*). That
;                 is, the semi-major axis should be calculated
;                 from Kepler's equation: a/rstar = G(M+m)P^2/(4pi^2rstar).
;                 Using G = 2942.71377d0 R_sun^3/(m_sun*day^2), masses
;                 in solar masses, period in days, and rstar in solar
;                 radii will give the appropriate input.
;   tperiastron - periastron passage time (BJD)
;   Period      - Period of orbit (days)
;
; OPTIONAL INPUTS:  
;   e           - eccentricity of orbit (0 if not specified)
;   omega       - argument of periastron of the star's orbit, in radians
;                 -- omega_* is typically quoted from RV
;                 -- required if e is specified
;                 -- assumed to be pi/2 if e not specified
;                 -- omega_* = omega_planet + !dpi
;   lonascnode  - The Longitude of the ascending node
;                 (radians). Assumed to be !dpi if not specified.
;
; OUTPUTS:
;    result     - the impact parameter as a function of BJD, in units
;                 of R_*.
;
; OPTIONAL OUTPUTS:
;    x,y,z - Arrays of the cartesian coordinates of the planet at
;            each BJD, in the units of R_*. 
;              +X is right
;              +Y is up
;              +Z is out of the page (primary transit)
;
; MODIFICATION HISTORY 
;  2009/05/01 -- Jason Eastman (Ohio State University)
;  2017/04/12 -- Add optional mass ratio keyword
;-

function exofast_getb, bjd, i=i, a=a, tperiastron=tperiastron, Period=P, $
                       e=e, omega=omega, x=x, y=y, z=z, $
                       lonascnode=lonascnode, q=q

if n_elements(q) eq 0 then q = !values.d_infinity

;; calculate the mean anomaly corresponding to each observed time
meananom = (2.d0*!dpi*(1.d0 + (bjd - Tperiastron)/P)) mod (2.d0*!dpi) 

;; if eccentricity is given, integrate the orbit
if keyword_set(e) then begin
    if n_elements(omega) eq 0 then message, $
      'ERROR: omega must be specified if e is specified'
    eccanom = exofast_keplereq(meananom, e)
    trueanom = 2d0*atan(sqrt((1d0 + e)/(1d0 - e))*tan(0.5d0*eccanom))
endif else begin
    trueanom = meananom
    e=0.d0
endelse

atot = a ;; a1 + a2

;; calculate the corresponding (x,y) coordinates of planet
r = atot*(1d0-e^2)/(1d0+e*cos(trueanom))

;; as seen from observer
x = -r*cos(trueanom + omega)
tmp = r*sin(trueanom + omega)
y =  -tmp*cos(i)
z =  tmp*sin(i)

;; Rotate by the Longitude of Ascending Node
;; For transits, it is not constrained, so we assume Omega=!dpi)
if n_elements(lonascnode) eq 1 then begin
    xold = x & yold = y
    x = -xold*cos(lonascnode) + yold*sin(lonascnode)
    y = -xold*sin(lonascnode) - yold*cos(lonascnode)
endif

b = sqrt(x^2 + y^2)

return, b

end
 
