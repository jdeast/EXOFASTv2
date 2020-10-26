;+
; NAME:
;   EOPDATA
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Read and interpolate tabulated earth orientation parameters
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   EOPDATA, JDUTC, PMX, PMY, UT1_UTC, DPSI, DEPS, $
;        /RESET, FILENAME=, ANGUNITS=, TBASE=
;
; DESCRIPTION:
;
;  The procedure EOPDATA reads, interpolates and returns Earth
;  orientation parameters used for precision earth-base astronomy
;  applications.
;
;  ** NOTE: The user is responsible for downloading and maintaining an
;  up-to-date file of earth orientation parameters from the
;  International Earth Rotation Service.  See below. **
;
;  This interface is somewhat provisional.  See OPEN QUESTIONS below.
;
;  The values returned are described below.  These descriptions are
;  taken from the Explanatory Supplement to IERS Bulletins A and B.
;
;    * PMX and PMY, the coordinates of the Celestial Ephemeris Pole
;      (CEP) relative to the earth-fixed International Reference Pole
;      (IRP).  The x-axis is in the direction of the IERS Reference
;      Meridian (IRM), the y-axis is in the direction 90 degrees West
;      longitude.  The time series of PMX and PMY is referred to as
;      "polar motion."
;
;      These are the coordinates of the earth rotation pole, as seen
;      in an *earth-fixed* coordinate system.  A station whose
;      coordinates are given in earth-fixed coordinates referred to
;      the ITRS can be transformed to the earth-fixed coordinates
;      referred to the true rotation pole of date using the following
;      matrix transformation:
;            
;             R_TRUE = RX(PMY) ## RY(PMX) ## R_ITRS
;
;      where the matrices RX and RY are defined below.
;
;    * UT1, the the rotation angle about the pole. It is related to
;      the Greenwich mean sidereal time (GMST) by a conventional
;      relationship (Aoki et al., 1982).  It gives access to the
;      direction of the International Reference Meridian IRM in the
;      ICRS, reckoned around the CEP axis. It is expressed as the
;      difference UT1-UTC.  Thus, the value of UT1 is computed as:
;
;          UT1 = UT1_UTC + UTC
;
;      where UTC is the UTC time, expressed in seconds.
;
;    * DPSI and DEPS are the offsets in longitude and obliquity of the
;      celestial pole with respect to its direction defined using the
;      conventional IAU precession/nutation theory.  An a priori
;      correction model is available in the IERS Conventions (1996),
;      (McCarthy, 1996).  The expressions to compute the nutation
;      angles are:
;
;          DEPS_TRUE = DEPS_1980 + DEPS    ;; Nutation in obliquity
;          DPSI_TRUE = DPSI_1980 + DPSI    ;; Nutation in longitude
;
;      where DPSI_1980 and DEPS_1980 are the nutation values
;      determined from the IAU 1980 Nutation Theory; and DPSI_TRUE and
;      DEPS_TRUE are the nutations to be used as arguments to further
;      precession and nutation computations.
;
;  For requested times which are between tabular values, a linear
;  interpolation is performed.  This is not exactly the correct
;  procedure, and can result in errors of +/- 0.1 mas in the earth
;  polar motion and 1 usec in UT1 (see McCarthy & Gambis 1997).
;
;
; DATA FILES and MAINTENANCE
;
;   The user is responsible for downloading and maintaining the earth
;   orientation parameters file as supplied by the IERS.  The format
;   of the files is the "Final" EOP data ASCII format.  They can be
;   downloaded here:
;
;      ftp://maia.usno.navy.mil/ser7/finals.all   ;; from May 1976-present
;      ftp://maia.usno.navy.mil/ser7/finals.data  ;; from Jan 1992-present
;   
;   The user must place this file in a known location, and in *at
;   least the first call*, this filename must be passed using the
;   FILENAME keyword.
;
;   EOPDATA will load the data once on the first call, and keep a
;   cached copy for subsequent calls.  On a daily basis the file will
;   be reloaded in case the quantities have been updated from the
;   server.  A reload of data can be forced using the RESET keyword.
;
; ROTATION MATRICES
;
;   The rotation matrices RX(T) and RY(T) mentioned above in relation
;   to polar motion are:
;
;      RX(T) =EQ= [[1,0,0], [0,cos(T),sin(T)], [0,-sin(T),cos(T)]]
;      RY(T) =EQ= [[cos(T),0,-sin(T)], [0,1,0], [sin(T),0,cos(T)]]
;      RZ(T) =EQ= [[cos(T),sin(T),0], [-sin(T),cos(T),0], [0,0,1]]
;
;   and are meant to be applied to a vector R as, RX(T) ## R.
;
;
; OPEN QUESTIONS
;
;   How will the transition to a new IERS EOP series be accomplished?
;   Using a keyword?
;
;   Should there be a quality flag?  The EOP file contains a
;   "predicted" flag, and also there are rows which contain no value
;   at all.  These should probably be flagged somehow.
;
;
; INPUTS:
;
;   JDUTC - a vector or scalar, the UTC time for which earth
;           orientation parameters are requested, expressed in Julian
;           Days.  The value of the keyword TBASE is added to this
;           quantity to arrive at the actual Julian date.
;
; OUTPUTS:
;
;   PMX, PMY - the earth-fixed angular coordinates of the celestial
;              ephemeris pole, measured in ANGUNITS units.
;
;   UT1_UTC - the value of UT1 - UTC, expressed in seconds.
;
;   DPSI, DEPS - the corrections to the IAU 1980 theory of Nutation,
;                for nutation in longitude and obliquity, expressed in
;                ANGUNITS units.
;
; KEYWORD PARAMETERS:
;
;   FILENAME - scalar string, on the first call, the name of the file
;              from which earth orientation parameters will be read.
;              Default value: (none)
;
;   TBASE - a fixed epoch time (Julian days) to be added to each value
;           of JDUTC.  Since subtraction of large numbers occurs with
;           TBASE first, the greatest precision is achieved when TBASE
;           is expressed as a nearby julian epoch, JDUTC is expressed
;           as a small offset from the fixed epoch.
;           Default: 0
;
;   ANGUNITS - scalar string, output units of angular parameters.
;              Possible values are 'ARCSEC' or 'RADIAN'.
;              Default value: 'RADIAN'
;
;   RESET - if set, forces EOP file to be re-read.
;
;
; EXAMPLE:
;
;
; SEE ALSO:
;
;   HPRNUTANG, TAI_UTC (Markwardt Library)
;   PRECESS, PREMAT, JPRECESS, BPRECESS (IDL Astronomy Library)
;
;
; REFERENCES:
;
;   Aoki, S., Guinot, B., Kaplan, G.H., Kinoshita, H., McCarthy, D.D.,
;     Seidelmann, P.K., 1982: Astron. Astrophys., 105, 359-361.
;
;   McCarthy, D. D. (ed.) 1996: IERS Conventions, IERS T.N. 21.
;     http://maia.usno.navy.mil/conventions.html
;
;   McCarthy, D. \& Gambis, D. 1997, "Interpolating the IERS Earth
;     Orientation Data," IERS Gazette No. 13, 
;     http://maia.usno.navy.mil/iers-gaz13
;     Instructions for high precision EOP data interpolation, not done
;     in this procedure.
;
;   Ray, J. & Gambis, D. 2001, "Explanatory Supplement to IERS
;     Bulletins A and B,"
;     http://hpiers.obspm.fr/iers/bul/bulb/explanatory.html
;
;     Explains meanings of earth orientation parameters used and
;     returned by this procedure.
;
;   Definition of Final EOP data format
;     ftp://maia.usno.navy.mil/ser7/readme.finals
;
; MODIFICATION HISTORY:
;   Written, 30 Jan 2002, CM
;   Documented, 14 Feb 2002, CM
;   Add default message, 01 Mar 2002, CM
;   More robust handling of input file, 10 Mar 2002, CM
;   Fix bug in interpolation of nutation and polar motion adjustments,
;     thanks to Tim Lister, 2014-10-09, CM
;
;  $Id: eopdata.pro,v 1.7 2014/10/20 21:36:16 cmarkwar Exp $
;
;-
; Copyright (C) 2002, 2014, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

pro eopdata_read, file, jd, pmx, pmy, ut1, dpsi, deps, status=status

  status = 0
  get_lun, unit
  openr, unit, file, error=err
  if err NE 0 then begin
      free_lun, unit
      message, 'ERROR: could not open '+file
  endif

  buffersize = 2048
  ngood = 0L

  jd = 0d & pmx = 0D & pmy = 0D & ut1 = 0D & dpsi = 0d & deps = 0d
  ss = strarr(buffersize)
  while NOT eof(unit) do begin
      ;; Read data from disk
      on_ioerror, TRIM
      readf, unit, ss

      ;; On the last pass we may get less than a full buffer's worth
      TRIM: 
      cc = (fstat(unit)).transfer_count
      if cc EQ 0 then goto, DONE
      if cc LT buffersize then ss = ss(0:cc-1)

      ;; Parse the parameters
      mjd1 = double(strcompress(strmid(ss,7,9),/remove_all))
      NEXT0: on_ioerror, NEXT1
      pmx1 = double(strcompress(strmid(ss,18,10),/remove_all))
      NEXT1: on_ioerror, NEXT2
      pmy1 = double(strcompress(strmid(ss,37,10),/remove_all))
      NEXT2: on_ioerror, NEXT3
      ut11 = double(strcompress(strmid(ss,58,11),/remove_all))
      NEXT3: on_ioerror, NEXT4
      dps1 = double(strcompress(strmid(ss,97,11),/remove_all))
      NEXT4: on_ioerror, NEXT5
      dep1 = double(strcompress(strmid(ss,116,11),/remove_all))
      NEXT5:

      jd = [jd, mjd1]
      pmx  = [pmx, pmx1] & pmy = [pmy, pmy1]
      ut1  = [ut1, ut11]
      dpsi = [dpsi, dps1]
      deps = [deps, dep1]

      ss(*) = ''
  endwhile

  DONE:
  free_lun, unit
  if n_elements(jd) EQ 1 then begin
      message, 'ERROR: could not read data from '+file
  endif

  AS2R = !dpi/180d/3600d ;; Arcsec to radians
  jd = jd(1:*)
  pmx = pmx(1:*) * AS2R  ;; Convert to radians
  pmy = pmy(1:*) * AS2R  ;; Convert to radians
  ut1 = ut1(1:*)         ;; Already in seconds
  dpsi = dpsi(1:*) * (0.001D * AS2R) ;; Convert from mas to radians
  deps = deps(1:*) * (0.001D * AS2R) ;; Convert from mas to radians
  
  status = 1
  return
end

pro eopdata, jdutc, pmx, pmy, ut1_utc, dpsi, deps, reset=reset, $
             filename=filename0, angunits=angunits0, tbase=tbase0

  common eopdata_table, mjd0, pmx0, pmy0, ut10, dpsi0, deps0, $
    ntable, mjdmin, mjdmax, leap0, $
    timestamp, oldfile

  if n_params() EQ 0 AND n_elements(filename0) EQ 0 then begin
      USAGE:
      message, 'USAGE:', /info
      message, 'EOPDATA, JDUTC, PMX, PMY, UT1_UTC, DPSI, DEPS, '+$
        '[FILENAME=filename, ANGUNITS=angunits, TBASE=tbase', /info
      message, "   ANGUNITS is one of 'ARCSEC' or 'RADIAN'", /info
      return
  endif

  if n_elements(mjd0) EQ 0 OR keyword_set(reset) OR $
    (n_elements(jdutc) EQ 0 AND n_elements(filename0) GT 0) then begin
      RELOAD_COMMON:
      forward_function get_xtecal

      ;; Find filename
      if n_elements(filename0) EQ 0 then begin
          filename = ''
          ;; First try: the old file
          if n_elements(oldfile) NE 0 then $
            filename = oldfile

          ;; Second try: use standard IDL Astronomy Library data directory
          if filename EQ '' then $
            filename = find_with_def('iers_final_a.dat','ASTRO_DATA')

          ;; Third try: Markwardt-specific
          if filename EQ '' then begin
              catch, catcherr
              if catcherr EQ 0 then $
                filename = get_xtecal() + 'clock/finals.data'
              catch, /cancel
          endif

          ;; Could not find it, so trigger a fatal error
          if filename EQ '' then $
            message, 'ERROR: Could not find EOP data file'
      endif else begin
          filename = strtrim(filename0(0),2)
      endelse
      
      eopdata_read, filename, mjd0, pmx0, pmy0, ut10, dpsi0, deps0, $
        status=status

      ;; Make a series of pseudo leap seconds so that we can
      ;; interpolate linearly below.
      wh = where(abs(ut10(1:*) - ut10) GT 0.8 AND ut10 NE 0, ct)
      leap0 = lonarr(n_elements(ut10))
      if ct GT 0 then begin
          wh = [wh, n_elements(ut10)-1] & ct = ct + 1
          for i = 0, ct-2 do $
            leap0(wh(i)+1:wh(i+1)) = (leap0(wh(i)) $
                                      - round(ut10(wh(i)+1)-ut10(wh(i))))
      endif

      if status NE 1 then begin
          message, 'ERROR: could not read EOP data from '+filename
      endif

      timestamp = systime(1)
      oldfile = filename

      if n_elements(jdutc) EQ 0 then return
  endif

  if n_params() EQ 0 then goto, USAGE

  if systime(1) - timestamp GT 86400d then goto, RELOAD_COMMON

  if n_elements(tbase0) EQ 0 then $
    tbase = 0d $
  else $
    tbase = double(tbase0)

  if n_elements(angunits0) EQ 0 then $
    angunits = 'RADIAN' $
  else $
    angunits = strtrim(strupcase(strcompress(angunits0(0))),2)
  
  ;; Convert from Julian days to MJD
  t = (jdutc(*) + (tbase - 2400000.5d))

  tmin = min(t, max=tmax)
  wh = where(mjd0 GE tmin-5 AND mjd0 LE tmax+5 AND ut10 NE 0, ct)
  if ct EQ 0 then begin
      OUT_OF_BOUNDS:
      message, 'ERROR: input time was out of bounds'
  endif

  mjd1 = mjd0(wh)
  ii = value_locate(mjd1, t)
  jj = wh(ii)

  ;; This is UT1 - UTC(TSTART), so it should be continuous.  Thus, we
  ;; can interpolate it.  The normal UT1-UTC series has discontinuities.
  ut1 = (ut10+leap0)(jj)

  ;; Linear interpolation
  dt = (t-mjd1(ii))/(mjd1(ii+1)-mjd1(ii))
  ut1_utc = ut10(jj) +  dt * (ut1(ii+1)-ut1(ii))

  ;; Interpolate DPSI and DEPS, the adjustments to the nutation angles
  dpsi = dpsi0(jj) + dt * (dpsi0(jj+1)-dpsi0(jj))
  deps = deps0(jj) + dt * (deps0(jj+1)-deps0(jj))

  ;; Polar motion parameters
  pmx = pmx0(jj) + dt * (pmx0(jj+1)-pmx0(jj))
  pmy = pmy0(jj) + dt * (pmy0(jj+1)-pmy0(jj))

  ;; Units conversions
  case angunits of 
      'ARCSEC': begin
          R2AS = 3600d*180d/!dpi ;; Radian to arcsec
          
          dpsi = dpsi * R2AS
          deps = deps * R2AS
          pmx  = pmx  * R2AS
          pmy  = pmy  * R2AS
      end
      
      'RADIAN': begin
          dummy = 1
      end

      else: begin
          message, 'ERROR: angular unit '+angunits+$
            ' was not recognized'
      end
  end

  return
end

