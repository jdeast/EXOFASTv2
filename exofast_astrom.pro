;+
;  BJDTDB   - The BJD_TDB of the observations. 
;  RA/DEC   - The ICRS coordinates of the target (scalar degrees), at EPOCH
;             (default J2000). See Fig 10 for required precision.
;  PMRA     - Proper motion in RA, in mas/year (scalar)
;             Must be in units of arc (as is typical for most
;             modern catalogs); i.e. PMRA = d(RA)/dt * cos(DEC)
;  PMDEC    - Proper motion in dec, in mas/year (scalar)
;             See Fig 13 for required precision.
;  PX       - Parallax of target, in mas (scalar) (See Fig 11)
;  RV       - Radial velocity of target (m/s)
;             ** only ~100 km/s precision required for cm/s precision
;             decade intervals for nearby disk star (see Fig 12)**
;  EPOCH    - The epoch of the coordinates in JD_TDB, default is
;             julday(1,1,2000,12) = 2451545d0 => J2000
;  OBSPOS   - A 3xNTIMES vector of Solar System Barycenter coordinates
;             of the observatory, in AU.
;  STARPOS  - A 3xNTIMES vector of Target Barycentric coordinates of
;             the star, in AU. Default is zero (no companions). This
;             will generally be calculated by EXOFAST_GETB2.
;  EARTH    - If set, will calculate the geocentric Earth position for
;             the observatory position
;  GAIA     - If set, will calculate the position of Gaia for the
;             observatory position (not yet supported)
;-
function exofast_astrom, bjdtdb, ra, dec, pmra, pmdec, px=px, rv=rv, epoch=epoch, tbase=tbase, obspos=obspos, earth=earth, gaia=gaia, au=au, starpos=starpos

ntimes = n_elements(bjdtdb)

if n_elements(au) eq 0 then au = 1.495978707d11 ;; meters
if n_elements(px) eq 0 then px = 0d0
if n_elements(rv) eq 0 then rv = 0d0
if n_elements(epoch) eq 0 then epoch = julday(1,1,2000,12) ;; J2000
if n_elements(tbase) eq 0 then tbase = 0d0
if n_elements(starpos) ne 3L*ntimes then starpos = 0d0 ;; no companion

;; get the position of the observatory at each BJD_TDB
if n_elements(obspos) ne 3L*ntimes then begin
   ephemfile = find_with_def('JPLEPH.405','ASTRO_DATA')
   mindate = min(bjdtdb+tbase,max=maxdate) - 1.d0 & maxdate+=1
   if keyword_set(gaia) then begin
      message, 'Gaia ephemeris not yet supported!!'
   endif else if keyword_set(earth) then begin
      ;; use the geocenter of the Earth
      JPLEPHREAD,ephemfile,pinfo,pdata,[mindate,maxdate], $
                 status=status, errmsg=errmsg
      JPLEPHINTERP,pinfo,pdata,bjdtdb,x_earth,y_earth,z_earth,$
                   /earth,posunits="AU", tbase=tbase     
      obspos = transpose([[x_earth],[y_earth],[z_earth]])
   endif else begin
      message, 'OBSPOS, /EARTH, or /GAIA must be specified'
   endelse
endif

;; Follows prescription of Lindegren, 2017
;; http://adsabs.harvard.edu/abs/2012A%26A...538A..78L
rarad = ra*!dpi/180d0 ;; radians
decrad = dec*!dpi/180d0 ;; radians
pmraradperday = pmra*!dpi/(365.25d0*3600d3*180d0) ;; radians/day
pmdecradperday = pmdec*!dpi/(365.25d0*3600d3*180d0) ;; radians/day
pxrad = px*!dpi/(3600d3*180d0) ;; radians
rvradperday = rv*pxrad/AU*86400d0 ;; radians/day
starposrad = starpos*pxrad ;; radians

;; eq 5
p = [            -sin(rarad),             cos(rarad),        0d0]#replicate(1d0,ntimes)
q = [-sin(decrad)*cos(rarad),-sin(decrad)*sin(rarad),cos(decrad)]#replicate(1d0,ntimes)
r = [ cos(decrad)*cos(rarad), cos(decrad)*sin(rarad),sin(decrad)]#replicate(1d0,ntimes)

u = r + ((bjdtdb##replicate(1d0,3))-epoch)*(p*pmraradperday + q*pmdecradperday + r*rvradperday) - pxrad*obspos ;+ starposrad2;; eq 4
u /= (sqrt(total(u^2,1))##replicate(1d0,3)) ;; normalize vector

;; convert back to RA/Dec
raout = atan(u[1,*],u[0,*])*180d0/!dpi 
decout = asin(u[2,*])*180d0/!dpi  

;; this can't be right in detail....
if n_elements(starpos) ne 1 then begin
   raout += starposrad[0,*]*180d0/!dpi*cos(dec*!dpi/180d0)
   decout += starposrad[1,*]*180d0/!dpi
endif

neg = where(raout lt 0d0)
if neg[0] ne -1 then raout[neg] += 360d0

return, [raout,decout]

end
