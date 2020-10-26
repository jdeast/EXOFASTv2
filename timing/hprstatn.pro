;+
; NAME:
;   HPRSTATN
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compute high precision earth station positions in inertial coordinates
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   HPRSTATN, JDTT, R_ITRF, R_ECI, V_ECI, $
;             [ /JPL, /USE_EOP, /NO_UT1, TBASE= ]
;
; DESCRIPTION:
;
;  The procedure HPRSTATN computes the coordinates and velocities of
;  an earth station in J2000 equatorial earth-centered inertial
;  coordinates (ECI).  This may be useful in any application where an
;  earthbound observatory is used to collect data on a non-terrestrial
;  phenomenon.
;
;  The user must have the following routines involved: HPRNUTANG;
;  Markwardt Quaternion Library; JPLEPHREAD and JPLEPHINTERP (if JPL
;  keyword is used); EOPDATA (if USE_EOP keyword is set); and TAI_UTC.
;  Also, the appropriate data files for TAI-UTC and Earth Orientation
;  Parameters must be installed.
;
;  The user must specify the position of the earth station in
;  earth-centered, earth-fixed, cartesian coordinates of the ITRF.
;  The Z-axis points to terrestrial north, the X-axis lies in the
;  terrestrial equator pointing towards the Greenwich meridian, and
;  the Y-axis forms the right handed coordinate system.  Any
;  positional units may be specified.
;
;  For the highest precision, the preferred method is to know the
;  coordinates from a direct IRTF reduction, i.e., via VLBI.  A
;  procedure for estimating cartesian coordinates based on geodetic
;  coordinates (lat, lon) referred to the WGS84 ellipsoid is given
;  below.
;
;  The values returned are the earth-centered inertial J2000
;  coordinates and velocities of the station.  All the effects of
;  earth rotation, precession, nutation, and polar motion are
;  potentially included.  The user has a choice of the kinds of
;  transformations that are included (see JPL, USE_EOP and NO_UT1
;  keywords).
;
;  The returned positional units are the same as the input units.  The
;  returned velocity units are (units of input) PER SECOND.
;
;  It is possible specify more than one time, or more than one station
;  position, or both.  If both more than one time and position are
;  specified, then there must be an equal quantity of both.
;
;
; CARTESIAN COORDINATES FROM GEODETIC COORDINATES
;
;  For example, for a station whose geodetic latitude LAT, longitude
;  LON (where positive is east) and elevation H specified relative to
;  the ellipsoid, the cartesian coordinates are:
;
;         R_ITRF = [ (A*C + H)*COS(LAT)*COS(LON), $
;                    (A*C + H)*COS(LAT)*SIN(LON), $
;                    (A*S + H)*SIN(LAT) ]
;
;  where for the WGS84 reference ellipsoid, the equatorial radius is
;  set to A = 6378.137 km, and the flattening factor F =
;  1/298.257223563, and
;
;      C = SQRT(COS(LAT)^2 + (1 - F)^2*SIN(LAT)^2)
;      S = (1 - F)^2 * C
;
;
; INPUTS:
;
;   JDTT - a scalar or N-vector, the TT epoch time(s) for which
;          station coordinates are to be computed.
;
;          For reference, JDTT = JDTAI + 32.184/86400d, where JDTAI is
;          the international atomic time measured in days.  The value
;          of the keyword TBASE is added to JDTT to arrive at the
;          actual Julian date.
;
;   R_ITRF - cartesian coordinates of earth station.  Either a
;            3-vector, or a 3xN array.  Units can be any positional
;            units.
;
;
; OUTPUTS:
;
;   R_ECI - upon output, the coordinates of the station at the given
;           time(s), referred to the earth-centered J2000 coordinate
;           system.  Either a 3-vector or 3xN array depending on the
;           input.  Units are the same as for R_ITRF.
;
;   V_ECI - upon output, the velocities of the station at the given
;           time(s), referred to the earth-centered J2000 coordinate
;           system.  Either a 3-vector or 3xN array depending on the
;           input.  Units are (units of R_ITRF) PER SECOND.
;
;
; KEYWORD PARAMETERS:
;
;   TBASE - a fixed epoch time (Julian days) to be added to each value
;           of JDTT.  Since subtraction of large numbers occurs with
;           TBASE first, the greatest precision is achieved when TBASE
;           is expressed as a nearby julian epoch, JDTT is expressed
;           as a small offset from the fixed epoch.  
;           Default: 0
;
;   JPL - if set, then the JPL ephemeris is used to compute nutation
;         angles.  Otherwise the series representation of HPRNUTANG is
;         used.
;
;   USE_EOP - if set, then use earth orientation parameters, returned
;             by EOPDATA, to further refine the station coordinates.
;             Otherwise, only precession and nutation are used.
;
;   NO_UT1 - if set, then do not use the UT1-UTC conversion.
;   NO_PRECESSION - disable precession calculation.
;   NO_NUTATION - disable nutation calculation.
;   NO_POLAR_MOTION - disable polar motion calculation.
;
;
; EXAMPLE:
;
;   ;; ITRF coordinates of Deep Space Network Antenna 63 (METERS)
;   R_DSN63 = [+4849092.647d, -0360180.569d, +4115109.113d]
;
;   ;; Time: 2000/01/01 01:30 TT
;   JDTT = jday(2000d,1,1) + 1.5/24
;
;   ;; Compute position of antenna in J2000 coordinate system using
;   ;; full Earth Orientation Parameters.
;   HPRSTATN, JDTT, R_DSN63, R_ECI, V_ECI, /USE_EOP
;
;
; SEE ALSO:
;
;   HPRNUTANG, TAI_UTC (Markwardt Library)
;   EOPDATA, JPLEPHREAD, JPLEPHINTERP
;
;
; REFERENCES:
;
;   McCarthy, D. D. (ed.) 1996: IERS Conventions, IERS T.N. 21.
;     http://maia.usno.navy.mil/conventions.html
;
;   Seidelmann, P.K. 1992, *Explanatory Supplement to the Astronomical
;     Almanac*, ISBN 0-935702-68-7
;
;
; MODIFICATION HISTORY:
;   Written, 6 May 2002, CM
;   Documented, 12 May 2002, CM
;   Corrected discussion of geodetic coordinates, 26 May 2002, CM
;   Add NO_POLAR_MOTION keyword; only compute V_ECI if the variable
;    is passed, 07 Mar 2007, CM
;   Save some memory by deleting variables that are no longer used,
;    19 Dec 2008, CM
;
;  $Id: hprstatn.pro,v 1.6 2008/12/19 13:38:58 craigm Exp $
;
;-
; Copyright (C) 2002, 2007, 2008, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

pro hprstatn, jdtt, r_iers, r_eci, v_eci, tbase=tbase0, $
              jpl=jpl, use_eop=use_eop, no_ut1=no_ut1, $
              no_precession=noprec, no_nutation=nonut, $
              no_polar_motion=nopm

  ;; Compute the earth-orientation angles
  hprnutang, jdtt, tbase=tbase0, zeta, theta, z, dpsi, deps, $
    mean_obliq=eps0, true_obliq=eps, $
    gas_time=gast, $
    polar_x=pmx, polar_y=pmy, use_eop=use_eop, no_ut1=no_ut1, $
    no_nutation=nonut, jpl=jpl

  ;; Standard value of rotational angular velocity of the earth (Aoki
  ;; et al 1982)
  omega = 7.2921151467d-5 ;; RADIAN/SEC

  if keyword_set(nopm) then begin
      pmx = 0
      pmy = 0
  endif

  ;; Rotation from terrestrial reference pole to the true celestial
  ;; pole of date, then application of earth spin using apparent
  ;; sidereal time.
  ;;              GAST   PMY   PMX
  qter = qteuler(['z',   'x',  'y'], $
                 +gast,  -pmy, -pmx)
  gast = 0 & pmx = 0 & pmy = 0

  ;; Position referred to true celestial pole of date
  r_tod = qtvrot(r_iers, temporary(qter))

  ;; Velocity referred to true celestial pole of date
  if arg_present(v_eci) then begin
      v_tod = r_tod*0
      v_tod(0,*) = -omega*r_tod(1,*)
      v_tod(1,*) = +omega*r_tod(0,*)
  endif

  if keyword_set(noprec) then begin
      zeta = 0d & theta = 0d & z = 0d 
  endif
  if keyword_set(nonut) then begin
      eps0 = 0d & dpsi = 0d & eps = 0d
  endif

  ;; Application of earth precession and nutation from true equator
  ;; and equinox of date, to mean equator and equinox of J2000
  ;; (celestial)
  ;;              Precession        Nutation 
  qcel = qteuler(['z','y','z',      'x','z','x'], $
                 -zeta, +theta, -z, +eps0, -dpsi, -eps)
  zeta = 0 & theta = 0 & z = 0 & eps0 = 0 & dpsi = 0 & eps = 0

  ;; Compute station position in Earth Centered Inertial coordinates
  DO_PRNUTROT:
  r_eci = qtvrot(temporary(r_tod), qcel)
  if arg_present(v_eci) then $
    v_eci = qtvrot(temporary(v_tod), temporary(qcel))

  return
end
