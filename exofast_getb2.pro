;+
; NAME:
;   EXOFAST_GETB2
;
; PURPOSE: 
;   This function returns the impact parameter as a function of BJD
;   given orbital elements of many planets, assuming keplerian orbits.
;   Optionally, the 3-space barycentric coordinates of the star
;   (x1,y1,z1) and planets (x2,y2,z2) or the stellar coordinates of the
;   planets (x0,y0,z0) can be returned.
; 
; CALLING SEQUENCE:
;   result = getb(jds, i=i,a=a,tperiastron=tperiastron,period=p,$
;                 e=e,omega=omega, x=x, y=y, z=z)
;
; INPUTS:
;   bjd         - Barycentric Julians dates for desired impact
;                 parameters (ideally in the target's
;                 barycentric frame; see BJD2TARGET)
;   i           - an NPLANETS array of planetary inclinations (radians)
;   a           - an NPLANETS array of semi-major axis, a = a1+a2 (in
;                 units of R_*). That is, the semi-major axis should be calculated
;                 from Kepler's equation: a/rstar = G(M+m)P^2/(4pi^2rstar).
;                 Using G = 2942.71377d0 R_sun^3/(m_sun*day^2), masses
;                 in solar masses, period in days, and rstar in solar
;                 radii will give the appropriate input.
;   tperiastron - an NPLANETS array of periastron passage times (same
;                 units as BJD)
;   Period      - an NPLANETS array of Period of orbit (days)
;
; OPTIONAL INPUTS:  
;   e           - an NPLANETS array of eccentricities (0 if not specified)
;   omega       - an NPLANETS array of arguments of periastron of the star's orbit, in radians
;                 -- omega_* is typically quoted from RV
;                 -- assumed to be pi/2 if e not specified
;                 -- omega_* = omega_planet + !dpi
;   lonascnode  - an NPLANETS array of Longitudes of the ascending node
;                 (radians). Assumed to be !dpi if not specified.
;   Q           - an NPLANETS array of mass ratios of the primary to secondary
;                 (M1/M2). If not specified, the mass ratio is assumed
;                 to be infinite. The returned coordinates are the
;                 motion of the companion with respect to the primary.
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
;  2018/07/06 -- Now properly handles 2D time array (for long cadence)
;-

function exofast_getb2, bjd, inc=inc, a=a, tperiastron=tperiastron, Period=P, $
                       e=e, omega=omega, x0=x0, y0=y0, z0=z0, x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2,$
                       lonascnode=lonascnode, q=q

sz = size(bjd)
if sz[0] le 1 then begin
   ntimes = sz[1]
   ninterp = 1
endif else if sz[0] eq 2 then begin
   ;; long cadence extra sampling
   ninterp = sz[2]
   ntimes = sz[1]
endif else message, 'Incompatible dimensions on BJD'

nplanets = n_elements(inc)

if n_elements(q) ne nplanets then q = dblarr(nplanets) + !values.d_infinity
if n_elements(e) ne nplanets then e = dblarr(nplanets)
if n_elements(omega) ne nplanets then omega = dblarr(nplanets) + !dpi/2d0

x1 = dblarr(ntimes,ninterp)
y1 = dblarr(ntimes,ninterp)
z1 = dblarr(ntimes,ninterp)

x2 = dblarr(nplanets,ntimes,ninterp)
y2 = dblarr(nplanets,ntimes,ninterp)
z2 = dblarr(nplanets,ntimes,ninterp)

x0 = dblarr(nplanets,ntimes,ninterp)
y0 = dblarr(nplanets,ntimes,ninterp)
z0 = dblarr(nplanets,ntimes,ninterp)

isinfinite = where(~finite(q),complement=isfinite)
a2 = a*0d0
a1 = a*0d0
if isinfinite[0] ne -1 then begin
   a2[isinfinite] = a[isinfinite] 
   a1[isinfinite] = 0d0
endif 

if isfinite[0] ne -1 then begin
   a2[isfinite] = a[isfinite]*q[isfinite]/(1d0+q[isfinite])
   a1[isfinite] = a2[isfinite]/q[isfinite]
endif

for i=0L, nplanets-1L do begin

   ;; calculate the mean anomaly corresponding to each observed time
   meananom = (2.d0*!dpi*(1.d0 + (bjd - Tperiastron[i])/P[i])) mod (2.d0*!dpi) 

   ;; if eccentricity is given, integrate the orbit
   if e[i] ne 0d0 then begin
      eccanom = exofast_keplereq(meananom, e[i])
      trueanom = 2d0*atan(sqrt((1d0 + e[i])/(1d0 - e[i]))*tan(0.5d0*eccanom))
   endif else begin
      trueanom = meananom
   endelse

   ;; calculate the corresponding (x,y) coordinates of planet in
   ;; barycentric coordinates
   r2 = a2[i]*(1d0-e[i]^2)/(1d0+e[i]*cos(trueanom))

   ;; as seen from observer
   x2[i,*,*] = -r2*cos(trueanom + omega[i])
   tmp = r2*sin(trueanom + omega[i])
   y2[i,*,*] =  -tmp*cos(inc[i])
   z2[i,*,*] =  tmp*sin(inc[i])

   ;; Rotate by the Longitude of Ascending Node
   ;; For transits, it is not constrained, so we assume Omega=!dpi)
   if n_elements(lonascnode) eq nplanets then begin
      xold = x2[i,*,*] & yold = y2[i,*,*]
      x2[i,*] = -xold*cos(lonascnode[i]) + yold*sin(lonascnode[i])
      y2[i,*] = -xold*sin(lonascnode[i]) - yold*cos(lonascnode[i])
   endif

   ;; calculate the star position in the barycentric frame
   r1 = a1[i]*(1d0-e[i]^2)/(1d0+e[i]*cos(trueanom))

   ;; rotate to observer's plane of reference
   x1tmp = -r1*cos(trueanom + omega[i] + !dpi)
   tmp = r1*sin(trueanom + omega[i] + !dpi)
   y1tmp = -tmp*cos(inc[i])
   z1 += tmp*sin(inc[i])

   ;; Rotate by the Longitude of Ascending Node
   ;; For transits, it is not constrained, so we assume Omega=!dpi)
   if n_elements(lonascnode) eq nplanets then begin
      x1 += -x1tmp*cos(lonascnode[i]) + y1tmp*sin(lonascnode[i])
      y1 += -x1tmp*sin(lonascnode[i]) - y1tmp*cos(lonascnode[i])
   endif else begin
      x1 += x1tmp
      y1 += y1tmp
   endelse

endfor

;; now convert to stellar frame (which is relevant for transits)
for i=0L, nplanets-1L do begin
   x0[i,*,*] = x2[i,*,*] - x1
   y0[i,*,*] = y2[i,*,*] - y1
   z0[i,*,*] = z2[i,*,*] - z1
endfor

b = sqrt(x0^2 + y0^2)

return, b

end
 
