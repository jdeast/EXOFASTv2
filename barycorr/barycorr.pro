;+
; NAME:
;   BARYCORR
;
; PURPOSE:
;   Converts a measured redshift (zmeas) to a barycenter-corrected
;   Radial velocity, as described by Wright & Eastman, 2014.
;
;   If you use this code, please cite  Wright & Eastman, 2014.
;
; CALLING SEQUENCE:
;   ;; minimal
;   rv = barycorr(jdutc, ra, dec, zmeas, obsname=) 
;
;   ;; for highest precision (cm/s)
;   ;; See Wright & Eastman Sec 7 for required precision on each input
;   rv = barycorr(jdutc, ra, dec, zmeas, pmra=, pmdec=, px=, rv=, r_itrf=)
;
; INPUTS:
;  JDUTC    - The full julian date in UTC, (scalar or vector), e.g.,
;             2450000d0. Must be double precision.
;  RA/DEC   - The coordinates of the target (scalar degrees), in EPOCH
;             (default J2000). See Fig 10 for required precision.
;  ZMEAS    - The measured redshift (e.g., the result of cross correlation
;             with template spectrum). Scalar or vector.
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
;  LAT      - Latitude of the observatory, in degrees (scalar)
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
;  RV       - Radial velocit of target, in m/s
;             ** only ~100 km/s precision required for cm/s precision
;             decade timescales for nearby disk star (see Fig 12)**
;
; OPTIONAL KEYWORDS:
;
;  HIP      - If set, use the Hipparcos epoch (overrides EPOCH keyword)
;  SKIP_EOP - If set, suppress warnings about
;             $ASTRO_DATA/iers_final_a.dat and skip the higher precision
;             values for the nutation and precession of the Earth.
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
;-
function barycorr, jdutc, ra, dec, zmeas, pmra=pmra, pmdec=pmdec, px=px, rv=rv, $
                    obsname=obsname,lat=lat,long=long,alt=alt,r_itrf=r_itrf, $
                    skip_eop=skip_eop, tbase=tbase, epoch=epoch,hip=hip,$
                    no_ut1=no_ut1, no_precession=no_precession, $
                    no_nutation=no_nutuation, no_polar_motion=no_polar_motion, nogravity=nogravity, galilean=galilean, nolighttravel=nolighttravel, relative=relative, tt=tt, expmeter=expmeter, jdutcmid=jdutcmid


c = 299792458d0 ;; speed of light, in m/s
zb = zbarycorr(jdutc, ra, dec, pmra=pmra, pmdec=pmdec, px=px, rv=rv, $
               obsname=obsname, lat=lat, long=long, alt=alt, r_itrf=r_itrf, $
               skip_eop=skip_eop, tbase=tbase, epoch=epoch,hip=hip,$
               no_ut1=no_ut1, no_precession=no_precession, $
               no_nutation=no_nutuation, no_polar_motion=no_polar_motion, nogravity=nogravity, galilean=galilean, nolighttravel=nolighttravel, relative=relative, tt=tt, expmeter=expmeter, jdutcmid=jdutcmid)

;print, zmeas*zb*c
;stop
if n_elements(zmeas) ne n_elements(zb) then zmeas = replicate(zmeas,n_elements(zb))
;forprint, zmeas[0:10], zb[0:10], ((1d0+zmeas)*(1d0+zb) - 1d0)[0:10],/t
;print, '****************************'

return, c*((1d0+zmeas)*(1d0+zb) - 1d0) ;; eq 10

end
