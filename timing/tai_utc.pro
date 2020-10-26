;+
; NAME:
;   TAI_UTC
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compute (TAI - UTC) time difference (i.e., leap seconds)
;
; MAJOR TOPICS:
;   Time
;
; CALLING SEQUENCE:
;   LEAP = TAI_UTC(JD_UTC)                  ;; or,
;   LEAP = TAI_UTC(JD_TAI, /INVERT)
;
; DESCRIPTION:
;
;   The function TAI_UTC computes the difference between International
;   Atomic Time (TAI) and Universal Coordinated Time (UTC), in
;   seconds.  
;
;   After 01 Jan 1972, the two time systems are synchronized, except
;   for a number of leap seconds added to account for the varying rate
;   of rotation of the earth.  While TAI is a continuous atomic time
;   system, UTC is a civil time system which may have discontinuities
;   where leap seconds are introduced.  This function computes the
;   differences between the two time systems.
;
;   The conversion from UTC to TAI is computed as:
;
;     JD_TAI = JD_UTC + TAI_UTC(JD_UTC)/86400
;
;   The inversion conversion, from TAI to UTC, is computed as:
;
;     JD_UTC = JD_TAI + TAI_UTC(JD_TAI, /INVERT)/86400
;
;   Here JD_UTC and JD_TAI are the UTC and TAI times measured in
;   Julian days respectively.
;
;   The introduction of leap seconds is not predictable, owing to the
;   non-linear processes that govern the rotation of the earth.  The
;   International Earth Rotation Service determines when leap seconds
;   will be introduced.  Thus, the user must download the history of
;   leap seconds.  This file can be downloaded at the following URL:
;
;     ftp://maia.usno.navy.mil/ser7/tai-utc.dat
;
;   NOTE - the leap second file must be kept up to date as new leap
;   seconds are introduced.  The file is kept internally in
;   memory, but is reloaded from disk at least once per day.
;
;   If the disk file is not available, then a copy of the file as
;   available from the USNO in 2009 is used, but a warning message is
;   printed.
;
;   The leap second data can be loaded in several ways:
;      1. The FILENAME keyword may specify the exact file name and path;
;      2. If FILENAME is not defined, or the empty string, then
;         the default location $ASTRO_DATA/tai-utc.dat is used;
;         (ASTRO_DATA is a system environment variable, used by
;          the IDL astronomy library to store auxiliary data files)
;      3. If neither #1 or #2 are available, then the internal table
;         is used.
;
;
; PARAMETERS: 
;
;   JD - time measured in Julian days.  The time being converted
;        *from*.
;
; RETURNS:
;
;   The number of seconds to be added to the input time, to arrive at
;   the desired time.
;
;
; KEYWORD PARAMETERS:
;
;   INVERT - if set, then convert from TAI to UTC.  If not set
;            (default), then convert from UTC to TAI.
;
;   FILENAME - a scalar string, indicating the file name containing
;              leap second data.  The data is only loaded once upon
;              the first call, and then with a frequency determined by
;              the RELOAD_EVERY keyword.  If FILENAME is not
;              specified or a blank string, then the leap second data
;              is found using the methods described above.
;              Default: not defined; i.e. TAI_UTC searches the default
;                       locations
;
;   RELOAD_EVERY - a scalar value, indicates how often the data should
;                  be reloaded for long-running tasks.  The value is
;                  expressed in days.  If the leap second data was
;                  loaded more than RELOAD_EVERY days ago, then it
;                  will be reloaded.  Note that a value of 0 will
;                  cause immediate re-load of data.
;                  Default: 1 (i.e. re-load every 1 day)
; 
; EXAMPLE:
;
;   For data stored in $ASTRO_DATA,
;     print, tai_utc(2451544.5d)  ;; Uses $ASTRO_DATA/tai-utc.dat
;         32.000000
;     
;   For the data stored in one's home directory,
;     filename = getenv('HOME')+'tai-utc.dat'
;     print, tai_utc(2451544.5d, filename=filename)
;         32.000000
;
;
; REFERENCES:
;
;   Definition of leap seconds.
;      http://tycho.usno.navy.mil/leapsec.html
;
;   File containing leap seconds.
;     ftp://maia.usno.navy.mil/ser7/tai-utc.dat
;
;
; SEE ALSO
;   TDB2TDT, SYSTIME, CALDAT, JULDAY
;   
; MODIFICATION HISTORY:
;   Written and Documented, CM, Dec 2001
;   Fixed array indexing errors when the requested time range falls in
;     the leap second period, and the input is an array; avoided use
;     of variable JDAY, which is a function clash for me, 02 Mar 2002,
;     CM
;   Added helpful usage message, CM, 15 Mar 2002
;   Made file handling more robust (instead of crashing), CM, 19 Jul 2005
;   Add 01 Jan 2006 leap second, CM, 03 Oct 2005
;   Add 01 Jan 2009 leap second, CM, 21 Jul 2008
;   Add documentation and the RELOAD_EVERY keyword, CM, 02 Dec 2009
;   New default file location is $ASTRO_DATA/tai-utc.dat, CM, 28 Dec 2009
;
;  $Id: tai_utc.pro,v 1.9 2009/12/03 02:12:34 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, 2005, 2008, 2009, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes, and to modify and use for personal or
; internal use, is granted.  All other rights are reserved.
;-

pro tai_utc_preload, strs, msg
  msg = 'last leap second 2009 JAN 1'
  strs = [ $
' 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S', $
' 1961 AUG  1 =JD 2437512.5  TAI-UTC=   1.3728180 S + (MJD - 37300.) X 0.001296 S', $
' 1962 JAN  1 =JD 2437665.5  TAI-UTC=   1.8458580 S + (MJD - 37665.) X 0.0011232S', $
' 1963 NOV  1 =JD 2438334.5  TAI-UTC=   1.9458580 S + (MJD - 37665.) X 0.0011232S', $
' 1964 JAN  1 =JD 2438395.5  TAI-UTC=   3.2401300 S + (MJD - 38761.) X 0.001296 S', $
' 1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S', $
' 1964 SEP  1 =JD 2438639.5  TAI-UTC=   3.4401300 S + (MJD - 38761.) X 0.001296 S', $
' 1965 JAN  1 =JD 2438761.5  TAI-UTC=   3.5401300 S + (MJD - 38761.) X 0.001296 S', $
' 1965 MAR  1 =JD 2438820.5  TAI-UTC=   3.6401300 S + (MJD - 38761.) X 0.001296 S', $
' 1965 JUL  1 =JD 2438942.5  TAI-UTC=   3.7401300 S + (MJD - 38761.) X 0.001296 S', $
' 1965 SEP  1 =JD 2439004.5  TAI-UTC=   3.8401300 S + (MJD - 38761.) X 0.001296 S', $
' 1966 JAN  1 =JD 2439126.5  TAI-UTC=   4.3131700 S + (MJD - 39126.) X 0.002592 S', $
' 1968 FEB  1 =JD 2439887.5  TAI-UTC=   4.2131700 S + (MJD - 39126.) X 0.002592 S', $
' 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S', $
' 1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S', $
' 1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S', $
' 1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S', $
' 1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S', $
' 1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S', $
' 1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S', $
' 1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S', $
' 1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S', $
' 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S', $
' 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S', $
' 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S', $
' 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S', $
' 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S', $
' 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S', $
' 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S', $
' 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S', $
' 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S', $
' 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S', $
' 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S', $
' 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S', $
' 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S', $
' 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S', $
' 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S', $
' 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S']
  return
end

function tai_utc, jd, reset=reset, invert=invert, filename=filename0, $
                  reload_every=reload_every0

  common tai_utc_common, taiutc, timestamp

  if n_params() EQ 0 OR n_elements(jd) EQ 0 then begin
      message, 'USAGE: ', /info
      message, ' TAI = UTC + TAI_UTC(JD_UTC)      ;; (or)', /info
      message, ' UTC = TAI + TAI_UTC(JD_TAI, /INV)', /info
      message, '  ;; (JD_UTC is Julian date referred to UTC)', /info
      message, '  ;; (JD_TAI is Julian date referred to TAI)', /info
      message, '', /info
      message, 'Other timescales (all units in seconds; JD=Julian date):', $
        /info
      message, ' TT  = TAI + 32.184d     ;; TAI to Terrestrial Time', /info
      message, ' TDT = TAI + 32.184d     ;; TAI to Terrestrial Dynamical Time',$
        /info
      message, ' TDB = TDT + TDB2TDT(JD) ;; TDT to Barycentric Dynamical Time',$
        /info
      message, '     ;; (JD referred to TDT or TDB)', /info
      message, ' TAI = UTC + TAI_UTC(JD) ;; UTC to TAI (JD referred to UTC)', $
        /info
      message, ' UTC = TAI + TAI_UTC(JD, /INV) ;; UTC to TAI (JD referred to TAI)', /info
      return, 0
  endif

  if n_elements(reload_every0) EQ 0 then reload_every = 1d $
  else reload_every = double(reload_every0(0))

  if n_elements(taiutc) EQ 0 OR keyword_set(reset) then begin
      ;; Markwardt-specific function
      forward_function get_xtecal

      RELOAD_COMMON:
      root = {tai_utc_struct, jday: 0D, taiutc0: 0D, taiutc1: 0D, mjdoff: 0D}

      ;; Check default locations
      ;;   1. FILENAME given explicitly
      ;;   2. FILENAME not given but $ASTRO_DATA/tai-utc.dat exists
      ;;   3. Markwardt get_xtecal() works
      sz = size(filename0)
      if sz(sz(0)+1) EQ 7 then efile = strtrim(filename0(0),2) $
      else efile = find_with_def('tai-utc.dat','ASTRO_DATA')

      if efile EQ '' then begin
          catch, catcherr
          if catcherr EQ 0 then efile = get_xtecal()+'clock/tai-utc.dat'
          catch, /cancel
      endif

      if efile EQ '' then $
        message, 'ERROR: could not find TAI-UTC offset file'

      ;; Read the file
      catch, catcherr
      if catcherr EQ 0 then begin
          get_lun, unit
;          print, 'FILENAME=',efile
          openr, unit, efile, error=err
          if err EQ 0 then begin
              ;; Read the data from the file
              strs = ['']
              nline = 0L
              while NOT eof(unit) do begin
                  line = ''
                  readf, unit, line
                  strs = [strs, line]
                  nline = nline + 1
              endwhile
              free_lun, unit
              if nline LE 0 then begin
                  message, 'WARNING: no data was found in '+efile, /info
                  return, -1D
              endif

              strs = strs(1:*)

          endif else begin
              free_lun, unit
              message, 'WARNING: could not open '+efile, /info
              catcherr = 1
          endelse
      endif
      catch, /cancel
      
      if keyword_set(catcherr) then begin
          ;; Finding the file failed - use built-in table
          tai_utc_preload, strs, msg
          nline = n_elements(strs)
          message, 'WARNING: using canned data ('+msg+')', /info
      endif

      PARSE_STRS:
      taiutc = replicate(root, nline)
      taiutc.jday    = double(strmid(strs, 16,10))
      taiutc.taiutc0 = double(strmid(strs, 36,12))
      taiutc.mjdoff  = double(strmid(strs, 59,7))
      taiutc.taiutc1 = double(strmid(strs, 69,10))
      
      timestamp = systime(1)
      reload_every = 1d
  endif

  if systime(1) - timestamp GT reload_every*86400d then goto, RELOAD_COMMON

  jj = value_locate(taiutc.jday, jd)
  ii = jj > 0
  leaps = taiutc(ii).taiutc0 + $
          taiutc(ii).taiutc1 * (jd - 2400000.5D - taiutc(ii).mjdoff)
  wh = where(jj LT 0, ct) 
  if ct GT 0 then leaps(wh) = 0

  if keyword_set(invert) then begin
      ;; Special case where the leap second pushes us over a boundary
      wh = where(jj GE 0 AND jd-leaps/86400d LT taiutc(ii).jday, ct)
      if ct GT 0 then begin
          ii = ii(wh) - 1
          leaps(wh) = taiutc(ii).taiutc0 + $
            taiutc(ii).taiutc1 * (jd(wh) - 2400000.5D - taiutc(ii).mjdoff)
      endif
      leaps = -leaps
  endif

  return, leaps
end
