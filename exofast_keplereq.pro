FUNCTION exofast_keplereq,m,ecc,thresh=thresh
;+
; NAME:
;    exofast_keplereq
; PURPOSE: 
;    Solve Kepler's Equation
; DESCRIPTION:
;    Solve Kepler's Equation. Method by S. Mikkola (1987) Celestial
;       Mechanics, 40 , 329-334. 
;    result from Mikkola then used as starting value for
;       Newton-Raphson iteration to extend the applicability of this
;       function to higher eccentricities
;
; CATEGORY:
;    Celestial Mechanics
; CALLING SEQUENCE:
;    eccanom=exofast_keplereq(m,ecc)
; INPUTS:
;    m    - Mean anomaly (radians; can be an array)
;    ecc  - Eccentricity
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;    thresh: stopping criterion for the Newton Raphson iteration; the
;            iteration stops once abs(E-Eold)<thresh
; OUTPUTS:
;    the function returns the eccentric anomaly
; KEYWORD OUTPUT PARAMETERS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;  2002/05/29 - Marc W. Buie, Lowell Observatory.  Ported from fortran routines
;    supplied by Larry Wasserman and Ted Bowell.
;    http://www.lowell.edu/users/buie/
;
;  2002-09-09 -- Joern Wilms, IAA Tuebingen, Astronomie.
;    use analytical values obtained for the low eccentricity case as
;    starting values for a Newton-Raphson method to allow high
;    eccentricity values as well
;
;  $Log: keplereq.pro,v $
;  Revision 1.4  2009/11/04 Jason Eastman (OSU)
;  Changed name from keplereq to exofast_keplereq.pro
;  Fix infinite loop when diff is near 2pi +/- epsilon.
;  Fix infinite loop when oldval and eccanom[i] differ by integer
;  multiples of 2pi.
;  Consistent implementation of THRESH
;  Made 1/6 and 1/24 take full advantage of IDL double precision
;  Vectorized for loop -- factors of several faster for (rare) cases
;  with many elements to be refined
;
;  Revision 1.3  2005/05/25 16:11:35  wilms
;  speed up: Newton Raphson is only done if necessary
;  (i.e., almost never)
;
;  Revision 1.2  2004/08/05 10:02:05  wilms
;  now also works for more than 32000 time values
;
;  Revision 1.1  2002/09/09 14:54:11  wilms
;  initial release into aitlib
;-

  ;; set default values
  IF n_elements(thresh) EQ 0 THEN thresh=1d-10

  IF (ecc LT 0. OR ecc GE 1.) THEN BEGIN 
      message,'Eccentricity must be 0<= ecc. < 1'
  ENDIF 

  ;; so we don't modify the inputs
  mx=m

  ;; Range reduction of m to -pi < m <= pi
  ;; ... m > pi
  zz=where(mx GT !dpi)
  IF zz[0] NE -1 THEN BEGIN 
      mx[zz]=mx[zz] MOD (2*!dpi)
      zz=where(mx GT !dpi)
      IF zz[0] NE -1 THEN mx[zz]=mx[zz]-2.0D0*!dpi
  ENDIF 

  ;; ... m < -pi
  zz=where(mx LE -!dpi)
  IF zz[0] NE -1 THEN BEGIN 
      mx[zz]=mx[zz] MOD (2*!dpi)
      zz=where(mx LE -!dpi)
      IF zz[0] ne -1 THEN mx[zz]=mx[zz]+2.0D0*!dpi
  ENDIF 

  ;; Bail out for circular orbits...
  IF (ecc EQ 0.) THEN return,mx

  ;; equation 9a
  aux   =  4.d0*ecc+0.5d0
  alpha = (1.d0-ecc)/aux
  beta=mx/(2.d0*aux)

  ;; equation 9b (except not really; is there an errata?)
  ;; the actual equation 9b is much much slower, but gives the same
  ;; answer (probably because more refinement necessary)
  aux=sqrt(beta*beta+alpha*alpha*alpha)
  z=beta+aux
  zz=where(z LE 0.0d0)
  if zz[0] ne -1 THEN z[zz]=beta[zz]-aux[zz]
  z = z^0.3333333333333333d0

  s0=z-alpha/z
  s1=s0-(0.078d0*s0*s0*s0*s0*s0)/(1.d0+ecc)
  e0=mx+ecc*(3.d0*s1-4.d0*s1*s1*s1)

  se0=sin(e0)
  ce0=cos(e0)

  f  = e0-ecc*se0-mx
  f1 = 1.d0-ecc*ce0
  f2 = ecc*se0
  f3 = ecc*ce0
  u1 = -f/f1
  u2 = -f/(f1+0.5d0*f2*u1)
  u3 = -f/(f1+0.5d0*f2*u2+0.16666666666666667d0*f3*u2*u2)
  u4 = -f/(f1+0.5d0*f2*u3+0.16666666666666667d0*f3*u3*u3-0.041666666666666667d0*f2*u3*u3*u3)
            
  eccanom=e0+u4

  zz = where(eccanom GE 2.0d0*!dpi)
  IF zz[0] NE -1 THEN eccanom[zz]=eccanom[zz]-2.0d0*!dpi
  zz = where(eccanom LT 0.0d0)
  IF zz[0] NE -1 THEN eccanom[zz]=eccanom[zz]+2.0d0*!dpi

  ;; Now get more precise solution using Newton Raphson method
  ;; for those times when the Kepler equation is not yet solved
  ;; to better than 1e-10
  ;; (modification J. Wilms)

  ndx=where(mx LT 0.)
  IF (ndx[0] NE -1 ) THEN mx[ndx]=mx[ndx]+2.*!dpi

  ;; calculate the differences
  diff = abs(eccanom-ecc*sin(eccanom) - mx)
  diff = diff < abs(diff - 2*!dpi) ;; 0 - epsilon = 2*pi - epsilon
  ndx = where(diff GT thresh)

  WHILE NDX[0] NE -1 DO BEGIN
      ;; E-e sinE-M
      fe=(eccanom[ndx]-ecc*sin(eccanom[ndx])-mx[ndx]) MOD (2*!dpi)
      ;; f' = 1-e*cosE
      fs=(1.d0-ecc*cos(eccanom[ndx])) MOD (2*!dpi)
      oldval=eccanom[ndx]
      eccanom[ndx]=(oldval-fe/fs)

      tmp = where(abs(oldval-eccanom[ndx]) GT thresh)
      if tmp[0] ne -1 then ndx = ndx[tmp] $
      else break
  ENDWHILE
  ;; range reduction
  toohigh = where(eccanom GE 2*!dpi)
  if toohigh[0] ne -1 then eccanom[toohigh] = eccanom[toohigh] MOD (2*!dpi)
  toolow = where(eccanom LT 0)
  if toolow[0] ne -1 then eccanom[toolow] = (eccanom[toolow] MOD (2*!dpi)) + 2*!dpi 

  return,eccanom
END 
