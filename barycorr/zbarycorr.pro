;+
; NAME:
;   ZBARYCORR
;
; PURPOSE:
;   Determine the barycentric redshift (z_B) for a given star, as
;   described by Wright & Eastman, 2014. All equation numbers are
;   referenced to that paper. This function is intended to be used in
;   conjunction with BARYCORR.
;
;   If you use this code, please cite  Wright & Eastman, 2014.
;
; CALLING SEQUENCE:
;
;   ;; minimal
;   zb = zbarycorr(jdutc, ra, dec, obsname=) 
;
;   ;; for highest precision (cm/s)
;   ;; See Wright & Eastman Sec 7 for required precision on each input
;   zb = zbarycorr(jdutc, ra, dec, pmra=, pmdec=, px=, rv=, r_itrf=)
;
; INPUTS:
;  JDUTC    - The full Julian date in UTC, (scalar or vector), e.g.,
;             2450000d0. Must be double precision.
;  RA/DEC   - The coordinates of the target (scalar degrees), at EPOCH
;             (default J2000). See Fig 10 for required precision.
;
; OPTIONAL INPUTS:
;
;  OBSNAME  - The name of the observatory, as an input to
;             OBSERVATORY.PRO to retreive the latitude, longitude, and
;             altitude of the observing station. Either this,
;             latitude, longitude, and altitude, or R_ITRF must be
;             specified. If set, lat, long, alt, and r_itrf are
;             ignored.
;             See Fig 14 for required precision.
;  LAT      - Latitude (geodetic) of the observatory, in degrees (scalar)
;  LONG     - East longitude of the observatory, in degrees (scalar)
;  ALT      - Altitude of the observatory, in meters (scalar)
;  R_ITRF   - Three element array containing the XYZ coordinates of
;             observatory relative to geocenter, in meters. Only used
;             if obsname, lat, long, and alt are not specified.
;  EPOCH    - The epoch of the coordinates in JD_TDB, default is
;             julday(1,1,2000,12) = 2451545d0 => J2000
;             Overridden by HIP keyword.
;  TBASE    - The baseline that has been subtracted from the input
;             JD_UTCs, for higher precision times.
;
;             **** THESE ARE REQUIRED FOR CM/S PRECISION ****
;  PMRA     - Proper motion in RA, in mas/year (scalar)
;             Must be in units of arc (as is typical for most
;             modern catalogs); i.e. PMRA = d(RA)/dt * cos(DEC)
;  PMDEC    - Proper motion in dec, in mas/year (scalar)
;             See Fig 13 for required precision.
;  PX       - Parallax of target, in mas (scalar) (See Fig 11)
;  RV       - Radial velocity of target, in m/s
;             ** only ~100 km/s precision required for cm/s precision
;             decade intervals for nearby disk star (see Fig 12)**
;
; OPTIONAL KEYWORDS:
;
;  HIP      - If set, use the Hipparcos epoch (overrides EPOCH keyword)
;  SKIP_EOP - If set, suppress warnings about
;             $ASTRO_DATA/iers_final_a.dat and skip the higher precision
;             values for the Nutation and precession of the Earth.
;             **Do not set for cm/s precision; see Fig 7**
;
;  *** Use of these keywords is not recommended ***
;  NO_UT1          - Passed to HPRSTATN; if set, disable the UT1-UTC conversion. 
;  NO_PRECESSION   - Passed to HPRSTATN; if set, disable precession calculation.
;  NO_NUTATION     - Passed to HPRSTATN; if set, disable nutation calculation.
;  NO_POLAR_MOTION - Passed to HPRSTATN; if set, disable polar motion calc.
;
; RESULT    - The barycentric redshift, z_B, as defined in Wright &
;             Eastman, 2014. To correct RVs to high precision, you
;             must use this and the measured redshift:
;             RV_true = c*((1d0+z_meas)*(1d0+z_B)-1d0)
;
;             Using only required parameters, this correction is good
;             to ~10 cm/s. With all optional parameters, the correction
;             is good to ~2 mm/s and is dominated by Earth rotation
;             effects (See Figs 2 & 3, Wright & Eastman, 2014).
;
; DEPENDENCIES: 
;    1) You must have an environment variable, $ASTRO_DATA, that
;    specifies the path containing:
;       a) JPLEPH.405 (http://idlastro.gsfc.nasa.gov/ftp/data/JPLEPH.405)
;       b) tai-utc.dat (automatically retrieved by JDUTCTOJDTDB)
;       c) iers_final_a.dat 
;          Add these lines in your crontab to keep it up to date:
; 5 1 * * * . ~/.bashrc ; wget -NP wget -NP $ASTRO_DATA ftp://maia.usno.navy.mil/ser7/finals.all 
; 6 1 * * * . ~/.bashrc ; cp $ASTRO_DATA/finals.all $ASTRO_DATA/iers_final_a.dat
;    2) EXOFAST library (http://astroutils.astronomy.ohio-state.edu/exofast/)
;    3) Goddard astro library (https://github.com/wlandsman/IDLAstro)
;    4) Markwardt library (http://www.physics.wisc.edu/~craigm/idl/)
;
; EXAMPLES:  
;
; MODIFICATION HISTORY
; 
;  2014/07 -- Public Release
;               Jason Eastman (LCOGT)
;               Jason Wright (Penn State)
;  2014/08 -- Minor clarifications of comments.
;  2015/04 -- *** NOT BACKWARD COMPATIBLE!! ***
;             ** Change to *East* longitude for consistency with other
;             EXOFAST codes! 
;             Only matters if you use the long keyword
;             obsname, r_itrf behavior unchanged.
;             *** NOT BACKWARD COMPATIBLE!! ***
;-

function zbarycorr, jdutc, ra, dec, pmra=pmra, pmdec=pmdec, px=px, rv=rv, $
                    obsname=obsname,lat=lat,long=long,alt=alt,r_itrf=r_itrf, $
                    skip_eop=skip_eop, tbase=tbase, epoch=epoch,hip=hip,$
                    no_ut1=no_ut1, no_precession=no_precession, $
                    no_nutation=no_nutuation, no_polar_motion=no_polar_motion, nogravity=nogravity, galilean=galilean, nolighttravel=nolighttravel, relative=relative, tt=tt, expmeter=expmeter, jdutcmid=jdutcmid

;; constants corresponding to GM for each of the solar system bodies
;; the product GM is known better than G & M individually
;; (http://ssd.jpl.nasa.gov/?constants), in order of JPLEPHINTERP
GMsun = 1.32712440018d20 ;; G*M_sun, in m^3/s^2 
X = 328900.56d0 ;; M_sun/(M_earth+M_moon)
Y = 81.30059d0  ;; M_earth/M_moon
GM = GMsun/[0d0,6023600d0,408523.71d0,(X*(Y+1d0))/Y,3098708d0,1047.3486d0,$
            3497.898d0,22902.98d0,19412.24d0,1.35d8,(X*(Y+1d0)),1d0]
c = 299792458d0 ;; speed of light in m/s
AU = 1.495978707d11 ;; AU in meters
year = 365.25d0*3600d0*24d0 ;; year in seconds
pctoau = 3600d0*180d0/!dpi ;; pc in AU
kmstoauyr = year/AU*1d3 ;; km/s in AU/yr

;; if TBASE not set, don't add anything to the JD_UTCs
if n_elements(tbase) eq 0 then tbase = 0d0

;; default to stationary in space (~10 cm/s error)
if n_elements(pmra) eq 0 then pmra = 0d0
if n_elements(pmdec) eq 0 then pmdec = 0d0
if n_elements(px) eq 0 then px = 0d0
if n_elements(rv) eq 0 then rv = 0d0
;if n_elements(nogravity) eq 0 then nogravity = -1
if n_elements(epoch) eq 0 then epoch = julday(1,1,2000,12,0,0)
if keyword_set(hip) then epoch = 2448348.56250d0

ntimes = n_elements(jdutc)
nexpmeter = n_elements(expmeter)
if nexpmeter ne 0 then begin
   if nexpmeter ne ntimes then $
      message, "ERROR: The number of exposure meter values must be equal to the number of input times"
endif

;; If SKIP_EOP not set, make sure iers_final_a.dat is correctly installed
if not keyword_set(SKIP_EOP) then begin
   USE_EOP=1
   eopname=file_search(getenv('ASTRO_DATA'),/mark_directory)+'iers_final_a.dat'
   fileinfo = file_info(eopname)
   if fileinfo.exists then begin
      fileage = systime(/julian,/utc) - julday(1,1,1970,0,0,fileinfo.mtime) 
      if fileage gt 1.5d0 then begin
         print, 'WARNING: $ASTRO_DATA/iers_final_a.dat is out of date (' + $
                  strtrim(fileage,2) + $
                  ' days old); see dependences. May cause errors up to 2 cm/s.'
      endif
   endif else begin
      print, 'WARNING: $ASTRO_DATA/iers_final_a.dat does not exist; '+$
             'see dependencies. Will cause errors up to 2 cm/s'
      USE_EOP=0
   endelse 
endif else USE_EOP=0

;; Convert to JD_TDB
jdtdb = jdutc2jdtdb(jdutc, jd_tt=jd_tt, tbase=tbase)

if keyword_set(tt) then jdtdb = jdutc + 67.184d0/86400

;; retrieve the lat/lon/alt by name
if keyword_set(obsname) then begin
   observatory, obsname, obs
   lat = obs.latitude
   long = -obs.longitude
   alt = obs.altitude
endif 

;; position of an observing station w.r.t geocenter
if n_elements(lat) eq 1 and n_elements(long) eq 1 and $
   n_elements(alt) eq 1 then begin
   LAT1 = lat*!dpi/180d0
   LON1 = long*!dpi/180d0
   F = 1.d0/298.257223563d0
   CC = 1.d0/SQRT(COS(LAT1)^2 + (1.d0 - F)^2*SIN(LAT1)^2)
   S = (1.d0 - F)^2 * CC
   A = 6378137d0               ; m
   R_ITRF = [(A*CC + alt)*COS(LAT1)*COS(LON1),$
             (A*CC + alt)*COS(LAT1)*SIN(LON1),$
             (A*S  + alt)*SIN(LAT1)] ; m
endif else if n_elements(r_itrf) ne 3 then $
   message, 'Must specify either OBSNAME, R_ITRF, or LAT, LONG, and ALT'

;; apply rotation, precession, nutation, and polar motion to Earth
HPRSTATN, jd_tt, R_ITRF, R_ECI, V_ECI, TBASE=TBASE,USE_EOP=USE_EOP,/jpl, $
          no_ut1=no_ut1, no_precession=no_precession, no_nutation=no_nutuation,$
          no_polar_motion=no_polar_motion

;; position and velocity of geocenter w.r.t barycenter
jplfile = find_with_def('JPLEPH.405','ASTRO_DATA')
if jplfile EQ '' then message,'ERROR - Cannot find JPL ephemeris file' 
JPLEPHREAD,jplfile,pinfo,pdata,[long(jdtdb[0])-1,long(jdtdb[ntimes-1])+1]+tbase
JPLEPHINTERP, pinfo, pdata, jdtdb, x,y,z,vx,vy,vz, /EARTH,/VELOCITY, $
              VELUNITS = 'KM/S',tbase=tbase

;; position and velocity of observing station w.r.t barycenter
robs = r_eci + transpose([[x],[y],[z]])*1d3  ;; m
vgeo = transpose([[vx],[vy],[vz]])*1d3       ;; m/s

;; relativistic addition of velocities 
vobs = (v_eci + vgeo)/(1d0 + v_eci*vgeo/c^2) ;; m/s
beta_earth = vobs/c

;; apply proper motion, parallax, RV correction to coordinates 
;; See Astronomical Almanac 2014, pg B72-B74, B68

; Unit vector to star at EPOCH (q in AA)
r0hat = [cos(ra*!dpi/180d0)*cos(dec*!dpi/180d0), $
         sin(ra*!dpi/180d0)*cos(dec*!dpi/180d0), $
         sin(dec*!dpi/180d0)]#replicate(1d0,ntimes) 

up = [0d0, 0d0, 1d0] ;; eq 14
east = crossp(up, r0hat[*,0])
east = east / sqrt(total(east^2)) ;; unit East Vector (eq 15)
north = crossp(r0hat[*,0], east)  ;; unit North Vector (eq 16)
mu = ((pmra*east+pmdec*north)/pctoau/1d3)#replicate(1d0,ntimes) ;;rad/yr (eq 17)

;; stellar position at each time
epoch0 = 2000d0 + (epoch-2451545d0)/365.25d0 
yearnow = 2000d0 + (jdtdb+tbase-2451545d0)/365.25d0  
T = (yearnow-epoch0)##replicate(1d0,3)         ;; years
vpi = rv/1d3 * kmstoauyr * (px/1d3/pctoau)     ;; rad/yr
vel = mu + vpi*r0hat                           ;; rad/yr (m in AA)
r = r0hat + vel*T                              ;; rad    (p1 in AA)
rhat = r/(replicate(1d0,3)#sqrt(total(r^2,1))) ;; unitless

;; include Earth's motion
if px gt 0d0 then begin
   rho = 1d3*rhat/px*pctoau - robs/au                   ;; AU
   rhohat = rho/(replicate(1d0,3)#sqrt(total(rho^2,1))) ;; unitless
   r0 = 1d3/px*pctoau*au                                ;; m
   beta_star = r0*mu/c/year + rv*r0hat/c                ;; m/s

   ;; Butkevich & Lindegren, 2014, eq D.11
   ;; http://arxiv.org/abs/1407.4664
   zlighttravel = rv*r0*total(mu[*,0]^2)/c^2*T[0,*]/year

endif else begin
   rhohat = rhat ;; unitless
   beta_star = 0d0
   zlighttravel = 0d0
endelse
if keyword_set(nolighttravel) then zlighttravel = 0d0

;; Calculate gravitational redshift due to Solar System bodies
;; only Sun (11) and Earth (3) are important at the cm/s level
calcndx = [11,3,5,6,10,2,7,8,4,1,9] ;; in order of maximum correction

sum = 0d0
zshapiro = 0d0
for i=0, n_elements(calcndx)-1 do begin
   if calcndx[i] eq 3 then xmag = sqrt(total(r_eci^2,1)) $
   else begin 
      JPLEPHINTERP,pinfo,pdata,jdtdb,x,y,z,objectname=calcndx[i],tbase=tbase
      x = robs - transpose([[x],[y],[z]])*1d3 ;; m
      xmag = sqrt(total(x^2,1)) ;; m
      xhat = x/(replicate(1d0,3)#xmag) ;; unitless
      ;; add the Shapiro delay for each body while we're at it (eq 27)
      zshapiro += -2d0*GM[calcndx[i]]/c^2/(xmag*(1d0 + total(xhat*rhohat,1)))*$
        total(beta_earth*(rhohat-total(xhat*rhohat,1)##replicate(1d0,3)*xhat),1)
   endelse
   sum += GM[calcndx[i]]/xmag ;; (m/s)^2
endfor

zgravity = 1d0/(1d0 + sum/c^2) - 1d0 ;; eq 2
if keyword_set(nogravity) then zgravity=0d0

;; determine the Barycentric RV correction (eq 25)
gamma_earth = 1d0/(sqrt(1d0-total(beta_earth^2,1)))
if keyword_set(relative) then gamma_earth = 1d0

zb = -zlighttravel -zshapiro - 1d0 + gamma_earth*$
     ((1d0+total(beta_earth*rhohat,1))*(1d0+total(beta_star*r0hat,1)))/$
     ((1d0+zgravity)*(1d0+total(beta_star*rhohat,1)))

if keyword_set(galilean) then begin
;   stop
   zb = total(beta_earth*rhat,1)
endif

if nexpmeter ne 0 then begin
   ;; it's easy to lose precision here 
   ;; normalize the exposure meter and subtract a zero point from JD_UTC
   expmeternorm = expmeter/mean(double(expmeter))
   jd0 = min(jdutc)

   zb = total(expmeternorm*zb)/total(expmeternorm)
   jdutcmid = total(expmeternorm*(jdutc-jd0))/total(expmeternorm) + jd0
endif


return, zb

end
