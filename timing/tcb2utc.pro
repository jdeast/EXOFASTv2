;+
; NAME:
;   BJDTDB2JDUTC
; PURPOSE:
;   Converts a JD in Barycentric Coordinate Time (JD_TCB) to a JD in
;   Universal Coordinate Time (JD_UTC) with an accuracy of ~1 us.
;
;   Since the ephemeris files on which the corrections are based use
;   Terestial Time (TT), which we're trying to compute, we must
;   iteratively call UTC2TCB to converge on the UTC that produces the
;   input TCB. Iterates ~5 times.
;
;   *** SEE UTC2TCB.PRO FOR IMPORTANT NOTES ON ACCURACY AND
;   CLARIFICATION ON THE DIFFERENCE BETWEEN BJD AND JD_TCB  ***
;
; INPUTS:
;   JD_TCB     - A scalar or array of JDs (in TBC). Must be double
;                precision.
;   RA/DEC     - A scalar specifying the Right Ascension and
;                Declinatoin of the object, in decimal degrees. J2000
;                is assumed unless /B1950 is set
;   
; OPTIONAL INPUTS:;
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
;                possible precision, TBASE = floor(jd) is recommended.
;
; OPTIONAL KEYWORDS:
;   B1950      - If set, the input coordinates are assumed to be in equinox 
;                B1950 coordinates. Default is the J2000 equinox.
;   TIME_DIFF  - To return the difference in seconds instead (either
;                to have the offset or for higher precision), use this
;                keyword. 
;                NOTE: For highest possible precision (~10 us),
;                TBC2UTC is called again with the determined UTC and
;                /TIME_DIFF keyword rather than returning
;                (UTC-TBC)*86400.d0. If speed is more important than
;                ~10 us precision, should calculate this yourself.
;                UTC = TCB + tcb2utc(jd,ra,dec,/time_diff)/86400.d0
;   UPDATE     - For accurate results, current empirical data must
;                be used. Set this keyword to spawn wget to download
;                the most recent data file (tai-utc.dat). This is
;                only necessary when leap seconds are added.
; DEPENDENCIES:
;   IDL astronomy library
;
;   UTC2TCB (Jason Eastman - OSU)
;
;   an environment variable "ASTRO_DATA", which specifies that path to
;   tai-utc.dat and JPLEPH.405, which can be found here:
;   ftp://maia.usno.navy.mil/ser7/tai-utc.dat
;   http://www.physics.wisc.edu/~craigm/idl/down/JPLEPH.405
;
;   Craig Markwardt's routines:
;   http://www.physics.wisc.edu/~craigm/idl/down/
;   TAI_UTC
;   (for earth position correction)
;   HPRSTATN HPRNUTANG QTEULER QTCOMPOSE QTMULT QTVROT
;
; REVISION HISTORY:
; 2009/12/1: Written by Jason Eastman (OSU)

function tcb2utc, tcb, ra, dec, lat, lon, elevation, obsname=obsname, $
                   B1950=b1950, tbase=tbase, TIME_DIFF=time_diff, UPDATE=update

;; --------------------------------------------------------
;; UPDATE THE DATA FILE (for leap seconds)
;; --------------------------------------------------------
if keyword_set(update) then $
  spawn, 'wget -qNP ' + getenv('ASTRO_DATA') + $
  ' ftp://maia.usno.navy.mil/ser7/tai-utc.dat' 

;; required for convergence (otherwise double precision isn't enough)
if n_elements(tbase) eq 0 then tbase = 0
base = floor(tcb)
tcbfloor = tcb - base
tbasetot = tbase + base

utc = tcbfloor - utc2tcb(tcbfloor,ra,dec,lat,lon,elevation,obsname=obsname, $
                    B1950=B1950, tbase=tbasetot, /time_diff)/86400.d0

;; iterative process to find the UTC that corresponds to the TCB
;; completes in ~5 iterations
repeat begin
    tcb_new =  utc2tcb(utc,ra,dec,lat,lon,elevation,obsname=obsname, $
                       B1950=B1950, tbase=tbasetot)
    diff = tcbfloor-tcb_new
    utc += diff
    
endrep until max(abs(diff)) eq 0 

;; return the difference in seconds, if desired (high precision)
if keyword_set(time_diff) then begin
    return, -utc2tcb(utc,ra,dec,lat,lon,elevation,obsname=obsname, $
                     B1950=B1950, tbase=tbasetot, /time_diff)  
endif

;; return UTC - TBASE (add my baseline back in)
return, utc + base


end    
    
