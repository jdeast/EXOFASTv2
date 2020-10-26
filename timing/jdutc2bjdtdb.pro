;+
; NAME:
;   JDUTC2BJDTDB
; PURPOSE:
;   Converts a coordinated universal time JD (JD_UTC) to a Barycentric
;   Dynamical Time BJD (BJD_TDB) with an accuracy of ~1 us. BJD_TDB is
;   often referred to by the less strictly defined "BJD" (see notes
;   below).
;
; CALLING SEQUENCE
;   BJD_TDB = JDUTC2BJDTDB(jd_utc, ra, dec [,lat, lon, elevation,
;                    obsname=obsname, tbase=tbase, /1950, /TIME_DIFF,
;                    /TT_IN])
;
; INPUTS:
;   JD_UTC     - A scalar or array of JDs (in UTC). Must be double
;                precision. 
;                NOTE: If the keyword TT_IN is set, must be an
;                array of JDs in TT.
;                JD_TT = JD_UTC + 32.184 + N
;                where N is the number of leap seconds since 1961.
;   RA/DEC     - A scalar specifying the Right Ascension and
;                Declinatoin of the object, in decimal degrees. J2000
;                is assumed unless /B1950 is set
;   
; OPTIONAL INPUTS:
;   LAT/LON    - The latitude and longitude of the observatory, in
;                decimal degrees. Longitude is west. If LAT/LON and
;                OBSNAME are not specified, the center of the
;                earth will be assumed and the resultant accuracy will
;                be ~20 us. If OBSNAME is set, these are ignored.
;   ELEVATION  - The elevation of the observatory, in meters. If
;                LAT/LON are not set, or if OBSNAME is set, this will
;                be ignored. Sea level is assumed if not given.
;   OBSNAME    - A string input to OBSERVATORY.PRO that specifies an
;                observatory. If set, LAT/LON/ELEVATION are ignored.
;   TBASE      - A scalar or n_elements(JD) vector that is the
;                baseline subtracted from the input JDs. For highest
;                precision, TBASE = floor(jd) is recommended
;
; OPTIONAL KEYWORDS:
;   B1950      - If set, the input coordinates are assumed to be in equinox 
;                B1950 coordinates. Default is the J2000 equinox.
;   TIME_DIFF  - To return the difference in seconds instead (either
;                to have the offset or for higher precision), use this
;                keyword. 
;                TDB = UTC + jdutc2bjdtdb(jd,ra,dec,/time_diff)/86400.d0
;   TT_IN      - Set this keyword to use JD_TT as an input
;                instead. This will skip the check for leap second
;                updates and assume you have done this already
; OUTPUTS:
;   BJD_TDB     - The light travel time corrected, Barycentric
;                 Dynamical Time (TDB) in BJD for each given UTC. If
;                 /TIME_DIFF is set, then returns the time difference
;                 in seconds.
;                TDB = jdutc2bjdtdb(jd_utc,ra,dec)
;                TDB = UTC + jdutc2bjdtdb(jd_utc,ra,dec,/time_diff)/86400.d0
; 
; DESCRIPTION:
;   TDB time is the current time standard (as of 2009). This routine
;   follows the procedure here:
;   http://lheawww.gsfc.nasa.gov/users/craigm/bary/
;   but ignores the dispersion correction. BJD_TDB is defined:
;   BJD_TDB = TOBS + GEOMETRIC + CLOCK + EINSTEIN - SHAPIRO
;
;   where
;
;     TOBS      - the observed JD in UTC of the event on Earth
;     GEOMETRIC - the light travel time from your position on earth to
;                 the barycenter of the solar system (~300 s)
;     CLOCK     - the correction from input time to TDB time
;     EINSTEIN  - the relativistic correction due to using the earth
;                 as your inertial frame (~1 ms)
;     SHAPIRO   - the time delay due to photon bending in the
;                 potential of the solar system (~1 us)
;
;   This correction has been verified to ~1 us accuracy and may be
;   accurate to ~30 ns. However, great care must be exercised to get
;   even ~1 minute accuracy.
;
;   *** PLEASE READ THE FOLLOWING NOTES ABOUT ACCURACY CAREFULLY ***
;
;   NOTES (in order of accuracy of result):
;
;   This routine transparently handles the conversion from UTC to TDB.
;   Using TT will produce a systematic offset equal to the number of
;   leapseconds required plus 32.184 seconds, unless the TT_IN keyword
;   is set (66.184 seconds in 2009). 
;
;   HELIO_JD.PRO, BARYCEN.PRO, and many routines that calculate the
;   position of astronomical objects assume the input "JD" is "JD_TT",
;   a critical assumption that is not obvious to users unfamiliar with
;   the complexities of precision timing. This has likely led to the
;   fact that the term "BJD" has been used in the literature to mean
;   both "BJD_UTC" and "BJD_TT" (similarly for "HJD"), which differ by
;   ~1 minute! BJD_UTC and HJD_UTC are not continous, uniformly
;   increasing timescales and should never be used in astronomy.
;
;   Leap seconds are added unpredicably every ~1 year. A current file
;   is vital for ~1 second accuracy. This program will wget an updated
;   file the first time it run2 after every January 1st and July 1st
;   (when leap seconds are added). Therefore wget and an internet
;   connection are required to use JD_UTC as an input. If wget fails,
;   you can manually retrieve it from
;   ftp://maia.usno.navy.mil/ser7/tai-utc.dat and put it in
;   $ASTRO_DATA. You'll also have to update
;   $ASTRO_DATA/exofast_lastupdate to have the current JD_UTC so this
;   program knows the leap seconds are up to date. To skip this step
;   altogether, input JD_TT and set the TT_IN keyword.
;
;   UTC and UT1 may differ by as much as 0.9s (the difference between
;   "delta t" and leap seconds + 32.184), and both may be called
;   UT. UTC is returned by NTP servers, and is usually the value
;   recorded in image headers, but this should be verified if UT is
;   specified.
;
;   For it to take into account your position on the earth, LAT and
;   LON, or OBSNAME must be specified. If none of these are set, the
;   assumed position is the center of the earth, and the results will
;   be systematically offset by 8-21.3 ms (observation bias). This
;   correction agrees with naive estimates (spherical earth at noon on
;   the winter solstice at -23.5 degrees latitude and 0 longitude) to
;   2 ms. Based on this and a handful of similar tests, this program
;   is believed to accurately implement Craig Marquardt's routines,
;   which are reported to be limited by the GPS coordinates of your
;   observatory.
;
;   ****************************************************************
;   The accuracy beyond this has not been thoroughly tested. While
;   routines such as BARYCEN and AXBARY claim an accuracy of better
;   than 1 us, they do not take into account the observer's position
;   on Earth, which is a 20 ms effect (4 orders of magnitude
;   worse!). For what it's worth, correcting to the geocenter, our
;   program agrees with BARYCEN to 60 ns, which agrees with FXBARY to
;   1 us. We believe ours is better than BARYCEN because it properly
;   indexes the JPL ephemeris with TDB, not TT. We also get consistent
;   results with GET_BJDTDB to 5 ns, though they are identical
;   programs save for the ephemeris interpolation. However, 
;   ***any accuracy beyond 1 ms should be independently verified.***
;   ****************************************************************
;
;   For 1 ms accuracy, double precision cannot hold the full
;   JD. TBASE=floor(jd) should be used.
;
;   ELEVATION or OBSNAME should be set to account for the light travel
;   time to your elevation. If not set, the elevation is assumed to be
;   at sea level and results will be accurate to ~10 us.
;
;   For ~100 ns accuracy, days cannot hold the precision. The
;   correction should be returned in seconds using the /TIME_DIFF
;   keyword and applied carefully.
;
;   For better than ~30 ns, the corrections are limited by the
;   imperfect knowledge of your location with respect to the center of
;   the earth (~10 m). A slight modification of this code to use ITRF
;   measurements of your position on earth (ie from VLBI) will be
;   required. Further verification of this routine at that level is
;   strongly recommended, and heroic measures to ensure double
;   precision is adequate may be required.
;
;   Also note the TT to TDB correction uses the 791-term Fairhead &
;   Bretagnon analytical approximation which has a maximum error of 23
;   ns and rms error of 14 ns in the time range 1980-2000, compared to
;   a numerical integration.
;
; DEPENDENCIES:
;   IDL astronomy library
;
;   wget and an internet connection (to update the leap second file)
;   alternatively, the $ASTRO_DATA/tai-utc.dat can be manually updated
;   and $ASTRO_DATA/exofast_lastupdate can be manually edited to
;   contain the JD_UTC of the manual update. Note this only happens
;   the first time it is run after every Jan 1st and Jul 1st.
;
;   an environment variable "ASTRO_DATA", which specifies that path to
;   JPLEPH.405, which can be found here:
;   http://www.physics.wisc.edu/~craigm/idl/down/JPLEPH.405
;   NOTE that tai-utc.dat and exofast_lastupdate will be placed in
;   $ASTRO_DATA too.
;   For manual updates, tai-utc.dat can be found here:
;   ftp://maia.usno.navy.mil/ser7/tai-utc.dat
;
;   Craig Markwardt's routines:
;   http://www.physics.wisc.edu/~craigm/idl/down/
;   TAI_UTC
;   (for earth position correction)
;   HPRSTATN HPRNUTANG QTEULER QTCOMPOSE QTMULT QTVROT
;
; REVISION HISTORY:
; 2010/04/08: Added TT_IN keyword
; 2009/12/10: Written by Jason Eastman (OSU)
; Based heavily on BARYCEN.PRO by goehler -- major differences are:
;   uses JD_UTC as input (not MJD_TT)
;   properly indexes JPL ephemeris with JD_TDB, not JD_TT (+/-60 ns effect)
;   automatically wgets of leap second file required to convert to TT
;   allows user settable TBASE for higher precision
;   uses your position on earth for higher precision (see notes)

function jdutc2bjdtdb, jd_utc, ra, dec, lat, lon, elevation, OBSNAME=obsname, $
                  TBASE=tbase, B1950=b1950, TIME_DIFF=time_diff,TT_IN=TT_IN
                  

ntimes = n_elements(jd_utc)
if n_elements(tbase) eq 0 then tbase = 0

;; convert to UTC to TDB 
;; (very close to Coordinate Time, the JPL ephemeris time)
jd_tdb = jdutc2jdtdb(jd_utc, TT_IN=TT_IN, tbase=tbase, $
                     jd_tt=jd_tt, clock_corr=clock_corr)

;; --------------------------------------------------------
;; read the ephemeris file with +/- 1 day margin
;; --------------------------------------------------------
ephemfile = find_with_def('JPLEPH.405','ASTRO_DATA')

IF NOT file_test(ephemfile) THEN $
  message,"Error: JPL Ephemeris file not found"

mindate = min(jd_tdb,max=maxdate) - 1.d0
maxdate += 1.d0
JPLEPHREAD,ephemfile,pinfo,pdata,[mindate,maxdate]+tbase, $
  status=status, errmsg=errmsg

IF status EQ 0 THEN message,"Ephemeris file reading failed: " + errmsg

;; --------------------------------------------------------
;; COMPUTE GEOMETRIC CORRECTION (~300 s)
;; --------------------------------------------------------
JPLEPHINTERP,pinfo,pdata,jd_tdb,x_earth,y_earth,z_earth,    $
  vx_earth, vy_earth,vz_earth,/earth,                      $
  posunits="KM", tbase=tbase,velunits='KM/S',/velocity
r_earth = transpose([[x_earth],[y_earth],[z_earth]])
rv_earth = transpose([[vx_earth],[vy_earth],[vz_earth]])

;; position on earth (~20 ms)
if n_elements(obsname) ne 0 then begin
    observatory, obsname, obs
    lat1 = obs.latitude*!dpi/180.d0
    lon1 = obs.longitude*!dpi/180.d0
    elevation = obs.altitude
endif else begin
    if n_elements(lat) eq 1 and n_elements(lon) eq 1 then begin
        lat1 = lat*!dpi/180.d0
        lon1 = lon*!dpi/180.d0
        if n_elements(elevation) eq 0 then elevation = 0
    endif
endelse

if n_elements(lat1) ne 0 then begin
    ;; calculate the cartesian coordinates of the observatory
    ;; described by the standard (geodetic) lat/lon
    F = 1.d0/298.257223563d0
    C = SQRT(COS(LAT1)^2 + (1 - F)^2*SIN(LAT1)^2)
    S = (1 - F)^2 * C
    H = elevation/1000.d0       ; km
    A = 6378.137                ; km
    R_ITRF = [(A*C + H)*COS(LAT1)*COS(LON1),$
              (A*C + H)*COS(LAT1)*SIN(LON1),$
              (A*S + H)*SIN(LAT1)]
    HPRSTATN, jd_tt, R_ITRF, R_ECI, V_ECI, TBASE=TBASE
    r_obs = r_earth + R_ECI
    einstein_corr = total(r_eci*rv_earth,1)/((pinfo.c / 1000.D0)^2)
endif else begin
    r_obs = r_earth
    einstein_corr = 0.d0
endelse

;; precess the object coordinates if desired
if keyword_set(B1950) then jprecess,ra,dec,ra1,dec1 else begin
    ra1 = ra
    dec1 = dec
endelse

;; convert coordinates to radians:
ra1 = ra1*!dpi/180.d0
dec1 = dec1*!dpi/180.d0

;; calculate the vector of the object
r_obj = [cos(dec1)*cos(ra1),$
         cos(dec1)*sin(ra1),$
         sin(dec1)]

;; scalar product of stellar object vector,earth vector, divided by c
;; in units of km/sec ( c is given in m/sec)
geo_corr = total(r_obs*(r_obj#replicate(1,ntimes)),1)/pinfo.c * 1000.d0

;; ------------------------------------------------------------
;; COMPUTE SHAPIRO CORRECTION (~20 us)
;; ------------------------------------------------------------
JPLEPHINTERP,pinfo,pdata,jd_tdb,x_sun,y_sun,z_sun,/sun,$
  posunits="KM", tbase=tbase
r_sun = transpose([[x_sun],[y_sun],[z_sun]])

;; distance of sun to observatory:
sun_dist = sqrt(total((r_sun - r_obs)^2,1))

;; cosine of unit vector sun->obs and unit vector obs -> object:
costh = total((r_obs-r_sun)*(r_obj#replicate(1,ntimes)),1)/sun_dist

;; apply shapiro correction. Sign in accordance with axbary.
;; refer  I.I. Shapiro, Phys. Rev. Lett. 13, 789 (1964)).
shapiro_corr = 2 *pinfo.msol * alog(1+costh)

;; ------------------------------------------------------------
;; SUMMARIZE CORRECTIONS
;; ------------------------------------------------------------
time = geo_corr + clock_corr + einstein_corr + shapiro_corr

;; return the desired value
if keyword_set(time_diff) then return, time
return, jd_utc + time/86400.d0

end
