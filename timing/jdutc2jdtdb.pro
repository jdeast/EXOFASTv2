;+
; NAME:
;   JDUTC2JDTDB
; PURPOSE:
;   Converts a coordinated universal time JD (JD_UTC) to a Barycentric
;   Dynamical Time JD (JD_TDB) with an accuracy of ~23 ns. Note this
;   does not do the geometric correction for the BJD.
;
; CALLING SEQUENCE
;   JD_TDB = JDUTC2JDTDB(jd_utc, tbase=tbase, /TT_IN, jd_tt=jd_tt, $
;                        clock_corr=clock_corr)
;
; INPUTS:
;   JD_UTC     - A scalar or array of JDs (in UTC). Must be double
;                precision. 
;                NOTE: If the keyword TT_IN is set, must be an
;                array of JDs in TT.
;                JD_TT = JD_UTC + 32.184 + N
;                where N is the number of leap seconds since 1961.
;   
; OPTIONAL INPUTS:
;   TBASE      - A scalar or n_elements(JD) vector that is the
;                baseline subtracted from the input JDs. For highest
;                precision, TBASE = floor(jd) is recommended
;   BIPMFILE   - The name of a file containing the BIPM-TAI
;                offsets from ftp://tai.bipm.org/TFG/TT%28BIPM%29/
;                Note that the previous year's file should be
;                downloaded, and the .ext file should be concatenated
;                at the end. As of 2010, this would be the TTBIPM.09
;                and TTBIPM09.ext files. This is only necessary for
;                30 us precision.
;
; OPTIONAL KEYWORDS:
;   TT_IN      - Set this keyword to use JD_TT as an input
;                instead. This will skip the check for leap second
;                updates and assume you have done this already
; OUTPUTS:
;   JD_TDB     - The JD in Barycentric Dynamical Time (JD_TDB) for
;                each given JD_UTC.
; OPTIONAL OUTPUTS:
;   JD_TT      - The JD in Terrestrial Time
;   CLOCK_CORR - The correction, in seconds from the input time to TDB.
;                jd_tdb = jd_utc + clock_corr/86400.d0
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
;   http://www.physics.wisc.edu/~craigm/idl/down/tai_utc.pro
;   http://www.physics.wisc.edu/~craigm/idl/down/tdb2tdt.pro

; REVISION HISTORY:
; 2010/04/08: Written by Jason Eastman (OSU)

function jdutc2jdtdb, jd_utc, TT_IN=TT_IN, tbase=tbase, $
                      jd_tt=jd_tt, clock_corr=clock_corr, bipmfile=bipmfile

if n_elements(tbase) eq 0 then tbase = 0

if not keyword_set(TT_IN) then begin
    taifile = find_with_def('tai-utc.dat','ASTRO_DATA')

    ;; --------------------------------------------------------
    ;; UPDATE THE DATA FILE (for leap seconds)
    ;; --------------------------------------------------------
    ;; the name of the file that contains the JD_UTC of last update
    updatefile = file_search(getenv('ASTRO_DATA'),/mark_directory) + $
      'exofast_lastupdate'

    ;; if the update file doesn't exist, create one that forces an update
    if not file_test(updatefile) then begin
        openw, updatelun, updatefile, /get_lun
        printf, updatelun, julday(1,1,1950)
        free_lun, updatelun
    endif

    ;; see when the last update was
    lastupdated = read_ascii(updatefile)
    now = systime(/julian, /utc)
    caldat, now, month, day, year
    caldat, lastupdated.field1, monthup, dayup, yearup 
    if (year gt yearup) or (month ge 6 and monthup lt 6) then begin
        ;; if Jan 1st or Jul 1st has passed without an update, update
        spawn, 'wget --timeout=5 -NP ' + getenv('ASTRO_DATA') + $
          ' ftp://maia.usno.navy.mil/ser7/tai-utc.dat', $
          output, count=nout, /stderr
        if strpos(output[nout-1],"Remote file no newer ") ne -1 or $
          strpos(output[nout-2],"saved") ne -1 then begin
            ;; if wget was successful, update the lastupdated file
            openw, updatelun, updatefile, /get_lun
            printf, updatelun, systime(/julian, /utc)
            free_lun, updatelun
        endif else begin
            ;; if wget was unsuccessful, print ERROR
            readcol, taifile, year, month, day, format='i,a,i',/silent
            nlines = n_elements(year)
            print, 'ERROR: Could not update leap second file.'
            print, 'Last leap second at ' + string(year[nlines-1],$
                     month[nlines-1],day[nlines-1], format='(i04,x,a3,x,i1)')
            print, 'If you do not have an internet connection, ' + $
              'manually update '+ $
              file_search(getenv('ASTRO_DATA'),/mark_directory) + taifile + $
              " and write today's JD_UTC to " + updatefile
            print, 'or input JD_TT and set the /TT_IN keyword'
            print, 'JD_TT = JD_UTC + 32.184 + N'
            print, 'where N is the number of leap seconds' 
            stop
        endelse
    endif

    ;; apply the BIPM-TAI correction if BIPM file supplied (30 us)
    if n_elements(bipmfile) ne 0 then begin
       readcol, bipmfile, mjd, junk, dt, format='d,d,d',/silent
       jd = mjd + 2400000.5d0 - tbase
       if max(jd) lt max(jd_utc) then begin
          nlines = n_elements(jd)
          caldat, jd[nlines-1] + tbase, month,day,year
          print, 'WARNING: BIPM offsets do not cover desired range.'
          message,'Type ".con" to extrapolate using most recent data from ' $
                  + string(year,month,day,format='(i04,"/",i02,"/",i02)')
       endif
       dtbipm = interpol(dt,jd,jd_utc,/quadratic)/1d6
    endif else dtbipm = 0
    
    ;; jd_utc to jd_tt offset
    clock_corr = tai_utc(jd_utc + tbase, filename=taifile) + 32.184d0 + dtbipm
        
endif else clock_corr = 0.d0 ;;input already in TT time

;; compute TT time
jd_tt = jd_utc + clock_corr/86400.d0

;; jd_utc to jd_tdb offset 
;; TDB2TDT is misnamed -- it's really TDT2TDB (TDT = TT)
;; it's also not vectorized
if n_elements(tbase) gt 1 then begin
   tmp = dblarr(n_elements(jd_tt))
   for i=0, n_elements(jd_tt)-1 do tmp[i] = TDB2TDT(jd_tt[i], tbase=tbase[i])
endif else tmp = TDB2TDT(jd_tt, tbase=tbase)
clock_corr += tmp

;; convert to TDB (same as JPL ephemeris time)
jd_tdb = jd_utc + clock_corr/86400.d0

return, jd_tdb

end
