;+
; NAME:
;   HPRNUTANG
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compute high precision earth precession, nutation and orientation angles
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   HPRNUTANG, JDTT, ZETA, THETA, Z, DPSI, DEPS, $
;              POLAR_X=PMX, POLAR_Y=PMY, JD_UT1=JD_UT1, /USE_EOPDATA, $
;              TBASE=, FIXED_EPOCH=, FIXED_BASE=, $
;              /JPL, /NO_UT1, $
;              MEAN_OBLIQUITY=EPS0, TRUE_OBLIQUITY=EPS, $
;              GMS_TIME=GMST, GAS_TIME=GAST, EQ_EQUIONOX=EQEQ
;              
;
; DESCRIPTION:
;
;  The procedure HPRNUTANG computes values of the earth
;  orientation-related angles, including precession and nutation,
;  which are used for high precision earth-based astronomy
;  applications.  
;
;  It is the goal of this procedure to provide all angles relevant in
;  determining the position of an earth station, as measured in an
;  earth-fixed coordinate system, and converting to space-fixed
;  coordinates.  This is useful in applications where observations by
;  a station in the earth-fixed frame are taken of an astrophysical
;  object which is in the non-rotating space-fixed frame.
;
;  This routine potentially depends on the following external
;  procedures, which also themselves depend on external data files:
;
;      EOPDATA - estimates Earth orientation parameters (if
;                USE_EOPDATA keyword is set), depends on earth
;                orientation data file.
;      TAI_UTC - computes time difference TAI - UTC (leap seconds),
;                depends on leap seconds file.
;
;  This interface is somewhat provisional.  See OPEN QUESTIONS below.
;
;  The user requests the quantities for a particular set of epoch
;  times, as measured in Julian days, in the system of Terrestrial
;  Dynamical Time ( = TDT = TT ).  
;
;  HPRNUTANG returns several quantities.  It is not possible to
;  describe each of these quantities in full detail in this
;  documentation.  The user is referred to the Explanatory Supplement
;  to the Astronomical Almanac (Sec 3.2) for more complete
;  descriptions.  The quantities are:
;
;    * ZETA, THETA, Z, which are euler angles representing the
;      precession of the mean celestial ephemeris pole with respect to
;      the space-fixed coordinate system defined by the FIXED epoch.
;      For a vector R_MEAN_OFDATE, whose space-fixed coordinates are
;      referred to the mean pole of date, the transformation to
;      space-fixed coordinates referred to the mean pole of the fixed
;      epoch is:
;
;        R_FIXED = qtvrot(R_MEAN_OFDATE, $
;                        qteuler(['z','y','z'], -zeta, +theta, -z))
;
;      By default the "fixed" epoch is J2000.0.  [ See below for
;      definitions of QTVROT and QTEULER. ]
;
;    * DPSI, DEPS, which are the angles representing the nutation in
;      longitude and obliquity of the true of-date celestial ephemeris
;      pole with respect to the mean pole of date.  For a vector
;      R_TRUE_OFDATE, whose space-fixed coordinates are referred to
;      the true pole of date, the transformation to space-fixed
;      coordinates referred to the mean pole of date is:
;
;        R_MEAN_OFDATE = qtvrot(R_TRUE_OFDATE, $
;                              qteuler(['x','z','x'], $
;                                      +eps0, -dpsi, -eps0-deps)
;
;      where EPS and EPS0 are defined below.
;
;    * EPS0, which is the mean obliquity of the ecliptic plane,
;      referred to the mean equator of date, at the requested epoch.
;      For a vector, R_ECL_OFDATE, whose space-fixed coordinates are
;      referred to the mean ecliptic and equinox of date, the
;      transformation to space-fixed coordinates referred to the mean
;      equator and equinox of date is:
;
;        R_MEAN_OFDATE = qtvrot(R_ECL_OFDATE, $
;                               qteuler(['x'], eps0)
;
;    * EPS, which is the true obliquity of the ecliptic plane,
;      referred to the mean equator of date, at the requested epoch.
;
;    * GMST, GAST, which are the mean and apparent Greenwich Sidereal
;      Times at the requested epoch.  For a vector R_TRUE_EARTHFIXED,
;      whose earth-fixed coordinates are referred to the true pole of
;      date, the transformation to space-fixed coordinates referred to
;      the true pole of date are:
;
;        R_TRUE_OFDATE = qtvrot(R_TRUE_EARTHFIXED, $
;                              qteuler(['z'], +gast))
;
;    * EQEQ, the equation of the equinoxes at the requested epoch.
;      This quantity may be more commonly known as the "precession of
;      the equinox."
;
;    * PMX, PMY, the coordinates of the celestial ephemeris pole as
;      measured in the earth-fixed coordinate system (set to zero if
;      the USE_EOPDATA keyword is not set).  For a vector
;      R_MEAN_EARTHFIXED, whose earth-fixed coordinates are referred
;      to the International Reference Pole, the transformation to
;      earth-fixed coordinates referred to the true pole of date are:
;
;        R_TRUE_EARTHFIXED = qtvrot(R_MEAN_EARTHFIXED, $
;                                  qteuler(['x','y'], -pmy, -pmx))
;
;      The vector R_MEAN_EARTHFIXED, could be for example, the
;      cartesian coordinates of a station on the earth, as determined
;      from its geodetic/geocentric latitude and longitude.
;
;    * JD_UT1, the UT1 time (expressed in Julian days) (set to UTC if
;      the USE_EOPDATA keyword is not set or if NO_UT1 is set).
;
;  Users may select different techniques to compute some of these
;  quantities.  See keywords JPL and USE_EOPDATA.
;
;
; OPEN QUESTIONS
;
;   How will the transition to a new IERS EOP series be accomplished?
;   Using a keyword?  How can users select different nutation series?
;   How can users select different fundamental arguments for the
;   planets.
;
;
; VERIFICATION
;
;   The precession and nutation quantities were compared against those
;   produced by the SLALIB telescope pointing library.  
;
;   For the epoch JD 2450449 (TT), the precession quantities of
;   HPRNUTANG agree numerically with SLALIB SLA_PREC to within 0.1
;   microarcseconds, and the nutation quantities agree SLALIB SLA_NUTC
;   to within 6 microarcseconds (and 54 microarcseconds in the mean
;   obliquity).  The GMST values agree with SLALIB SLA_GMSTA to better
;   than 1 nanosecond.  Of course this says nothing about the accuracy
;   of the IAU 1976/1980 precession and nutation models with respect
;   to the true precession and nutations.
;
;   The precession and nutation quantities computed in this procedure
;   -- ZETA, THETA, Z, DPSI and DEPS -- were also used to compute the
;   space-fixed coordinates of the Goldstone DSS-63 deep space network
;   tracking station.  These values were compared against values
;   produced by JPL Horizons ephemeris generator.  Agreement was found
;   at the 60 cm level.  Accuracy at that level is probably limited by
;   the JPL DE406 earth ephemeris used by Horizons.
;
;   Polar motion values were estimated at the same epoch using
;   EOPDATA, and applied to three orthogonal unit vectors.  The above
;   quaternion transformation produces the same coordinate values,
;   when compared against SLALIB_POLMO.
;
;
; QTEULER and QTVROT
;
;   The functions QTEULER and QTVROT are functions from the Markwardt
;   quaternion library.  QTEULER composes a chain of Euler-like
;   rotations into a single quaternion.  QTVROT applies a quaternion
;   rotation to a 3-vector.
;
;   The user need not use these functions.  Any function which
;   constructs a set of Euler-like rotations, and then applies them to
;   3-vectors will work fine.
;
;
; INPUTS:
;
;   JDTT - a vector or scalar, the TT epoch time(s) for which high
;          precision values are to be computed.
;
;          For reference, JDTT = JDTAI + 32.184/86400d, where JDTAI is
;          the international atomic time measured in days.  The value
;          of the keyword TBASE is added to JDTT to arrive at the
;          actual Julian date.
;
; OUTPUTS:
;
;   ZETA, THETA, Z - Euler angles of precession of the mean celestial
;                    ephemeris pole, expressed in ANGUNITS units.
;
;   DPSI, DEPS - the nutation angles in longitude and obliquity of the
;                true pole with respect to the mean pole, expressed in
;                ANGUNITS units.  By default the values are based on
;                the IAU 1980 theory of nutation.  The user can select
;                JPL to interpolate the JPL nutation ephemerides.
;
;                When USE_EOPDATA is set, the nutation angles are
;                augmented by the offset correction terms supplied in
;                the EOP file.
;
; KEYWORD PARAMETERS:
;
;   TBASE - scalar or vector, a fixed epoch time (Julian days) to be
;           added to each value of JDTT.  Since subtraction of large
;           numbers occurs with TBASE first, the greatest precision is
;           achieved when TBASE is expressed as a nearby julian epoch,
;           JDTT is expressed as a small offset from the fixed epoch.
;           Default: 0
;
;   FIXED_EPOCH - a scalar or vector number, the fixed epoch (in TT
;                 Julian Days) against which the precession angles of
;                 the mean pole are referred.  
;                 Default: JD 2451545.0 TT ( = J2000.0 )
;
;   FIXED_BASE - scalar or vector, a fixed epoch time to be added to
;                FIXED_EPOCH, in much the same way that TBASE is added
;                to JDTT.  Default: 0
;
;   POLAR_X, POLAR_Y - upon return, the quantities PMX and PMY, in
;                      ANGUNITS units, if USE_EOPDATA is set.  If
;                      USE_EOPDATA is not set then zero is returned
;                      for both PMX and PMY.
;
;   JD_UT1 - upon return, the time in the UT1 system at the requested
;            epoch, if the USE_EOPDATA keyword is set.  If the
;            USE_EOPDATA keyword is not set, or if NO_UT1 is set, then
;            the time in UTC is returned (which is guaranteed to be
;            within +/- 0.9 seconds of UT1).
;
;   MEAN_OBLIQUITY - upon return, the quantity EPS0, in ANGUNITS
;                    units.
;
;   TRUE_OBLIQUITY - upon return, the quantity EPS, in ANGUNITS units.
;
;   GMS_TIME - upon return, the quantity GMST in radians.
;
;   GAS_TIME - upon return, the quantity GAST in radians.
;
;   EQ_EQUINOX - upon return, the quantity EQEQ in ANGUNITS units.
;
;   ANGUNITS - scalar string, output units of angular parameters.
;              Possible values are 'ARCSEC' or 'RADIAN'.
;              Default value: 'RADIAN'
;
;   JPL - if set, then use the JPL nutation ephemeris to determine the
;         nutation angle quantities.
;
;   NO_UT1 - if set, then do not compute UT1, but use UTC instead.
;
;   USE_EOPDATA - if set, use the EOPDATA procedure to determine earth
;                 orientation parameters at the requested epoch.
;                 These include polar motion values, corrections to
;                 the 1980 IAU nutation theory, and the UT1
;                 correction.
;
;
; EXAMPLE:
;
;   Need an example converting topocentric to/from J2000.0
;
;   Need an example converting station position earth-fixed
;   coordinates to/from space-fixed coordinates.
;
;
;
; SEE ALSO:
;
;   HPRNUTANG, TAI_UTC (Markwardt Library)
;   PRECESS, NUTATE, PREMAT, JPRECESS, BPRECESS (IDL Astronomy Library)
;
;
; REFERENCES:
;
;   Aoki, S., Guinot, B., Kaplan, G.H., Kinoshita, H., McCarthy, D.D.,
;     Seidelmann, P.K., 1982: Astron. Astrophys., 105, 359-361.
;
;   HORIZONS, JPL Web-based ephemeris calculator (Ephemeris DE406)
;      http://ssd.jpl.nasa.gov/horizons.html
;
;   McCarthy, D. D. (ed.) 1996: IERS Conventions, IERS T.N. 21.
;     http://maia.usno.navy.mil/conventions.html
;
;   Seidelmann, P.K. 1992, *Explanatory Supplement to the Astronomical
;     Almanac*, ISBN 0-935702-68-7
;
;
; MODIFICATION HISTORY:
;   Written, 30 Jan 2002, CM
;   Documented, 15 Feb 2002, CM
;   Added docs about ecliptic; added default of 'RADIAN' to code; 01
;     Mar 2002, CM
;   Corrected equation of equinoxes (had DPSI*COS(EPS0), when it
;     should be DPSI*COS(EPS)), 01 Mar 2002, CM
;   Added default message, 04 Mar 2002, CM
;   Added more logic to detect JPL ephemeris file, 17 Mar 2002, CM
;   Corrected discussion of geodetic coordinates, 26 May 2002, CM
;   Documentation tweaks, 05 Jan 2004, CM
;   Some modifications to conserve memory, 22 Dec 2008, CM
;   Allow TBASE/FBASE to be a vector, 01 Jan 2009, CM
;
;  $Id: hprnutang.pro,v 1.16 2009/01/02 17:44:40 craigm Exp $
;
;-
; Copyright (C) 2002, 2004, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

pro hprnutang_init_iau1980, argfacts, psiamps, epsamps

; IAU 1980 Nutation Model
;   * taken from http://hpiers.obspm.fr/eop-pc/models/nut_iau1980
iau_1980_nut_strs = [ $
;  # Fortran format : (5I3,8X,6F10.1) 
;  #  Delaunay Arg. period(days)  
;  # lm ls F  D Om            psisin     t sin    eps cos     t cos
;   0    .    10   .    20        30        40        50         60     
   '  0  0  0  0  1 -6798.4 -171996.0    -174.2   92025.0       8.9     ', $
   '  0  0  2 -2  2   182.6  -13187.0      -1.6    5736.0      -3.1     ', $
   '  0  0  2  0  2    13.7   -2274.0      -0.2     977.0      -0.5     ', $
   '  0  0  0  0  2 -3399.2    2062.0       0.2    -895.0       0.5     ', $
   '  0 -1  0  0  0  -365.3   -1426.0       3.4      54.0      -0.1     ', $
   '  1  0  0  0  0    27.6     712.0       0.1      -7.0       0.0     ', $
   '  0  1  2 -2  2   121.7    -517.0       1.2     224.0      -0.6     ', $
   '  0  0  2  0  1    13.6    -386.0      -0.4     200.0       0.0     ', $
   '  1  0  2  0  2     9.1    -301.0       0.0     129.0      -0.1     ', $
   '  0 -1  2 -2  2   365.2     217.0      -0.5     -95.0       0.3     ', $
   ' -1  0  0  2  0    31.8     158.0       0.0      -1.0       0.0     ', $
   '  0  0  2 -2  1   177.8     129.0       0.1     -70.0       0.0     ', $
   ' -1  0  2  0  2    27.1     123.0       0.0     -53.0       0.0     ', $
   '  1  0  0  0  1    27.7      63.0       0.1     -33.0       0.0     ', $
   '  0  0  0  2  0    14.8      63.0       0.0      -2.0       0.0     ', $
   ' -1  0  2  2  2     9.6     -59.0       0.0      26.0       0.0     ', $
   ' -1  0  0  0  1   -27.4     -58.0      -0.1      32.0       0.0     ', $
   '  1  0  2  0  1     9.1     -51.0       0.0      27.0       0.0     ', $
   ' -2  0  0  2  0  -205.9     -48.0       0.0       1.0       0.0     ', $
   ' -2  0  2  0  1  1305.5      46.0       0.0     -24.0       0.0     ', $
   '  0  0  2  2  2     7.1     -38.0       0.0      16.0       0.0     ', $
   '  2  0  2  0  2     6.9     -31.0       0.0      13.0       0.0     ', $
   '  2  0  0  0  0    13.8      29.0       0.0      -1.0       0.0     ', $
   '  1  0  2 -2  2    23.9      29.0       0.0     -12.0       0.0     ', $
   '  0  0  2  0  0    13.6      26.0       0.0      -1.0       0.0     ', $
   '  0  0  2 -2  0   173.3     -22.0       0.0       0.0       0.0     ', $
   ' -1  0  2  0  1    27.0      21.0       0.0     -10.0       0.0     ', $
   '  0  2  0  0  0   182.6      17.0      -0.1       0.0       0.0     ', $
   '  0  2  2 -2  2    91.3     -16.0       0.1       7.0       0.0     ', $
   ' -1  0  0  2  1    32.0      16.0       0.0      -8.0       0.0     ', $
   '  0  1  0  0  1   386.0     -15.0       0.0       9.0       0.0     ', $
   '  1  0  0 -2  1   -31.7     -13.0       0.0       7.0       0.0     ', $
   '  0 -1  0  0  1  -346.6     -12.0       0.0       6.0       0.0     ', $
   '  2  0 -2  0  0 -1095.2      11.0       0.0       0.0       0.0     ', $
   ' -1  0  2  2  1     9.5     -10.0       0.0       5.0       0.0     '  ]
iau_1980_nut_strs = [ iau_1980_nut_strs, $
   '  1  0  2  2  2     5.6      -8.0       0.0       3.0       0.0     ', $
   '  0 -1  2  0  2    14.2      -7.0       0.0       3.0       0.0     ', $
   '  0  0  2  2  1     7.1      -7.0       0.0       3.0       0.0     ', $
   '  1  1  0 -2  0   -34.8      -7.0       0.0       0.0       0.0     ', $
   '  0  1  2  0  2    13.2       7.0       0.0      -3.0       0.0     ', $
   ' -2  0  0  2  1  -199.8      -6.0       0.0       3.0       0.0     ', $
   '  0  0  0  2  1    14.8      -6.0       0.0       3.0       0.0     ', $
   '  2  0  2 -2  2    12.8       6.0       0.0      -3.0       0.0     ', $
   '  1  0  0  2  0     9.6       6.0       0.0       0.0       0.0     ', $
   '  1  0  2 -2  1    23.9       6.0       0.0      -3.0       0.0     ', $
   '  0  0  0 -2  1   -14.7      -5.0       0.0       3.0       0.0     ', $
   '  0 -1  2 -2  1   346.6      -5.0       0.0       3.0       0.0     ', $
   '  2  0  2  0  1     6.9      -5.0       0.0       3.0       0.0     ', $
   '  1 -1  0  0  0    29.8       5.0       0.0       0.0       0.0     ', $
   '  1  0  0 -1  0   411.8      -4.0       0.0       0.0       0.0     ', $
   '  0  0  0  1  0    29.5      -4.0       0.0       0.0       0.0     ', $
   '  0  1  0 -2  0   -15.4      -4.0       0.0       0.0       0.0     ', $
   '  1  0 -2  0  0   -26.9       4.0       0.0       0.0       0.0     ', $
   '  2  0  0 -2  1   212.3       4.0       0.0      -2.0       0.0     ', $
   '  0  1  2 -2  1   119.6       4.0       0.0      -2.0       0.0     ', $
   '  1  1  0  0  0    25.6      -3.0       0.0       0.0       0.0     ', $
   '  1 -1  0 -1  0 -3232.9      -3.0       0.0       0.0       0.0     ', $
   ' -1 -1  2  2  2     9.8      -3.0       0.0       1.0       0.0     ', $
   '  0 -1  2  2  2     7.2      -3.0       0.0       1.0       0.0     ', $
   '  1 -1  2  0  2     9.4      -3.0       0.0       1.0       0.0     ', $
   '  3  0  2  0  2     5.5      -3.0       0.0       1.0       0.0     ', $
   ' -2  0  2  0  2  1615.7      -3.0       0.0       1.0       0.0     ', $
   '  1  0  2  0  0     9.1       3.0       0.0       0.0       0.0     ', $
   ' -1  0  2  4  2     5.8      -2.0       0.0       1.0       0.0     ', $
   '  1  0  0  0  2    27.8      -2.0       0.0       1.0       0.0     ', $
   ' -1  0  2 -2  1   -32.6      -2.0       0.0       1.0       0.0     ', $
   '  0 -2  2 -2  1  6786.3      -2.0       0.0       1.0       0.0     ', $
   ' -2  0  0  0  1   -13.7      -2.0       0.0       1.0       0.0     ', $
   '  2  0  0  0  1    13.8       2.0       0.0      -1.0       0.0     ', $
   '  3  0  0  0  0     9.2       2.0       0.0       0.0       0.0     ', $
   '  1  1  2  0  2     8.9       2.0       0.0      -1.0       0.0     ', $
   '  0  0  2  1  2     9.3       2.0       0.0      -1.0       0.0     ', $
   '  1  0  0  2  1     9.6      -1.0       0.0       0.0       0.0     ', $
   '  1  0  2  2  1     5.6      -1.0       0.0       1.0       0.0     ', $
   '  1  1  0 -2  1   -34.7      -1.0       0.0       0.0       0.0     ', $
   '  0  1  0  2  0    14.2      -1.0       0.0       0.0       0.0     ', $
   '  0  1  2 -2  0   117.5      -1.0       0.0       0.0       0.0     ', $
   '  0  1 -2  2  0  -329.8      -1.0       0.0       0.0       0.0     ', $
   '  1  0 -2  2  0    23.8      -1.0       0.0       0.0       0.0     ', $
   '  1  0 -2 -2  0    -9.5      -1.0       0.0       0.0       0.0     ', $
   '  1  0  2 -2  0    32.8      -1.0       0.0       0.0       0.0     ', $
   '  1  0  0 -4  0   -10.1      -1.0       0.0       0.0       0.0     ', $
   '  2  0  0 -4  0   -15.9      -1.0       0.0       0.0       0.0     ', $
   '  0  0  2  4  2     4.8      -1.0       0.0       0.0       0.0     ', $
   '  0  0  2 -1  2    25.4      -1.0       0.0       0.0       0.0     ', $
   ' -2  0  2  4  2     7.3      -1.0       0.0       1.0       0.0     ', $
   '  2  0  2  2  2     4.7      -1.0       0.0       0.0       0.0     ', $
   '  0 -1  2  0  1    14.2      -1.0       0.0       0.0       0.0     ', $
   '  0  0 -2  0  1   -13.6      -1.0       0.0       0.0       0.0     ', $
   '  0  0  4 -2  2    12.7       1.0       0.0       0.0       0.0     ', $
   '  0  1  0  0  2   409.2       1.0       0.0       0.0       0.0     ', $
   '  1  1  2 -2  2    22.5       1.0       0.0      -1.0       0.0     ', $
   '  3  0  2 -2  2     8.7       1.0       0.0       0.0       0.0     ', $
   ' -2  0  2  2  2    14.6       1.0       0.0      -1.0       0.0     ', $
   ' -1  0  0  0  2   -27.3       1.0       0.0      -1.0       0.0     ', $
   '  0  0 -2  2  1  -169.0       1.0       0.0       0.0       0.0     ', $
   '  0  1  2  0  1    13.1       1.0       0.0       0.0       0.0     ', $
   ' -1  0  4  0  2     9.1       1.0       0.0       0.0       0.0     ', $
   '  2  1  0 -2  0   131.7       1.0       0.0       0.0       0.0     ', $
   '  2  0  0  2  0     7.1       1.0       0.0       0.0       0.0     ', $
   '  2  0  2 -2  1    12.8       1.0       0.0      -1.0       0.0     ', $
   '  2  0 -2  0  1  -943.2       1.0       0.0       0.0       0.0     ', $
   '  1 -1  0 -2  0   -29.3       1.0       0.0       0.0       0.0     ', $
   ' -1  0  0  1  1  -388.3       1.0       0.0       0.0       0.0     ', $
   ' -1 -1  0  2  1    35.0       1.0       0.0       0.0       0.0     ', $
   '  0  1  0  1  0    27.3       1.0       0.0       0.0       0.0     ' ]

  n80 = n_elements(iau_1980_nut_strs)

  argfacts = fltarr(5,n80)
  str = strmid(iau_1980_nut_strs,0,16)
  reads, str, argfacts
  argfacts = transpose(argfacts)

  psiamps = dblarr(2,n80)
  str = strmid(iau_1980_nut_strs,23,21)
  reads, str, psiamps
  psiamps = transpose(psiamps)

  epsamps = fltarr(2,n80)
  str = strmid(iau_1980_nut_strs,43,20)
  reads, str, epsamps
  epsamps = transpose(epsamps)

  return
end

pro hprnutang_init_iau1980_args, args
  ;; c1 = mean anomaly of Moon
  ;; c2 = mean anomaly of Sun
  ;; c3 = mean longitude of the Moon minus the mean longitude of Moon's node
  ;; c4 = mean elongation of Moon from Sun
  ;; c5 = mean longitude of ascending node of the Moon
  ;; c6 = mean anomaly of Mercury
  ;; c7 = mean anomaly of Venus
  ;; c8 = mean anomaly of Earth
  ;; c9 = mean anomaly of Mars
  ;; ca = mean anomaly of Jupiter
  ;; cb = mean anomaly of Saturn

  c1 = [134.96298139d*3600d, 1717915922.6330d,  31.310d,   0.064d]
  c2 = [357.52772333d*3600d,  129596581.2240d,  -0.577d,  -0.012d]
  c3 = [ 93.27191028d*3600d, 1739527263.1370d, -13.257d,   0.011d]
  c4 = [297.85036306d*3600d, 1602961601.3280d,  -6.891d,   0.019d]
  c5 = [125.04452222d*3600d,   -6962890.5390d,   7.455d,   0.008d]
  c6 = [252.3d       *3600d,     149472.7d,      0d,       0d    ]
  c7 = [179.9d       *3600d,      58517.8d,      0d,       0d    ]
  c8 = [ 98.4d       *3600d,      35999.4d,      0d,       0d    ]
  c9 = [353.3d       *3600d,      19140.3d,      0d,       0d    ]
  ca = [ 32.3        *3600d,       3034.9d,      0d,       0d    ]
  cb = [ 48.0        *3600d,       1222.1d,      0d,       0d    ]

  args = [[c1],[c2],[c3],[c4],[c5],[c6],[c7],[c8],[c9],[ca],[cb]]
  args = args * !dpi / 180d / 3600d

  return
end

pro hprnutang_init_iau1996_args, args
  ;; c1 = mean anomaly of Moon
  ;; c2 = mean anomaly of Sun
  ;; c3 = mean longitude of the Moon minus the mean longitude of Moon's node
  ;; c4 = mean elongation of Moon from Sun
  ;; c5 = mean longitude of ascending node of the Moon
  ;; c6 = mean anomaly of Mercury
  ;; c7 = mean anomaly of Venus
  ;; c8 = mean anomaly of Earth
  ;; c9 = mean anomaly of Mars
  ;; ca = mean anomaly of Jupiter
  ;; cb = mean anomaly of Saturn
  ;; cc = accumulated general precession
  
  c1= [134.96340251d*3600d, 1717915923.2178d,  31.8792d,  5.1635d-2, 2.4470d-4]
  c2= [357.52910918d*3600d,  129596581.0481d,  -0.5532d,  1.36d-4,  -1.149d-5]
  c3= [ 93.27209062d*3600d, 1739527262.8478d, -12.7512d, -1.037d-3,  4.17d-6]
  c4= [297.85019547d*3600d, 1602961601.2090d,  -6.3706d,  6.593d-3, -3.169d-5]
  c5= [125.04455501d*3600d,   -6962890.2665d,   7.4722d,  7.702d-3, -5.939d-5]
  c6= [  0d,                     0d,                 0d,  0d,        0d      ]
  c7= [181.979800853d*3600d, 58517.8156748d*3600d,   0d,  0d,        0d      ]
  c8= [100.466448494d*3600d, 35999.3728521d*3600d,   0d,  0d,        0d      ]
  c9= [355.433274605d*3600d, 19140.299314d *3600d,   0d,  0d,        0d      ]
  ca= [ 34.351483900d*3600d,  3034.90567464d*3600d,  0d,  0d,        0d      ]
  cb= [ 50.0774713998d*3600d, 1222.11379404d*3600d,  0d,  0d,        0d      ]
  cc= [  0d,                  1.39697137214d*3600d,  3.086d-4, 0d,  0d       ]

  args = [[c1],[c2],[c3],[c4],[c5],[c6],[c7],[c8],[c9],[ca],[cb],[cc]]
  args = args * !dpi / 180d / 3600d

  return
end

pro hprnutang, jdtt, zeta, theta, z, dpsi, deps, jpl=jpl, $
               tbase=tbase0, polar_x=pmx, polar_y=pmy, $
               fixed_epoch=fepoch0, fixed_base=fbase0, $
               jd_ut1=jdut1, mean_obliquity=eps0, true_obliquity=eps, $
               gms_time=gmst, gas_time=gast, eq_equinox=eqeq, $
               use_eopdata=useeop, no_ut1=no_ut1, no_nutation=no_nut1

  common hprnutang_iau80_coeffs, arg80, psi80, eps80
  common hprnutang_iau80_args,   farg80
  common hprnutang_iau96_args,   farg96

  if n_params() EQ 0 then begin
      message, 'USAGE:', /info
      message, 'HPRNUTANG, JDTT, ZETA, THETA, Z, DPSI, DEPS, '+ $
        '[ TBASE=, /JPL, /USE_EOPDATA, /NO_UT1, FIXED_EPOCH=, FIXED_BASE=, '+ $
        'POLAR_X=, POLAR_Y=, JD_UT1=, MEAN_OBLIQUITY=, TRUE_OBLIQUITY=, '+$
        'GMS_TIME=, GAS_TIME=, EQ_EQUINOX= ]', /info
      return
  endif

  if n_elements(arg80) EQ 0 then begin
      hprnutang_init_iau1980, arg80, psi80, eps80 ;; IAU 1980 Nutation theory
      hprnutang_init_iau1980_args, farg80 ;; Fund. args of 1980 theory
      hprnutang_init_iau1996_args, farg96 ;; Fund. args of 1996 theory
  endif

; if keyword_set(arg96) then farg = farg96 else farg = farg80
  farg = farg80

  ;; Default angular units
  if n_elements(angunits0) EQ 0 then $
    angunits = 'RADIAN' $
  else $
    angunits = strtrim(strupcase(strcompress(angunits0(0))),2)

  ;; Default time bases
  if n_elements(tbase0) EQ 0 then tbase = 0d $
  else                            tbase = double(tbase0)

  if n_elements(fbase0) EQ 0 then fbase = 0d $
  else                            fbase = double(fbase0)
  
  ;; "Fixed" epoch, which is the epoch of equinox that coordinates are
  ;; precessed *to*, default 2000
  if n_elements(tfixed0) EQ 0 then fepoch = 2451545.0D - fbase $
  else                             fepoch = double(fepoch0(*))

  ;; Form epoch of date in centuries from J2000.0
  t = (jdtt(*) + (tbase - 2451545.0d))/36525d

  ;; Angular conversion factors
  TWOPI = 2d*!dpi
  AS2R  = !dpi/3600d/180d         ;; 1   arcsec to radians
  MAS2R = !dpi*0.0001d/3600d/180d ;; 0.1 milliarcsec to radians

  ;; Interpolate the JPL ephemerides of nutations if requested
  if keyword_set(jpl) then begin
      ;; Markwardt-specific function
      forward_function get_xtecal

      sz = size(jpl)
      if sz(sz(0)+1) EQ 7 then efile = strtrim(jpl(0),2) $
      else efile = find_with_def('JPLEPH.405','ASTRO_DATA')

      if efile EQ '' then begin
          catch, catcherr
          if catcherr EQ 0 then efile = get_xtecal()+'clock/JPLEPH.200'
          catch, /cancel
          if efile EQ '' then $
            message, 'ERROR: could not find JPL ephemeris'
      endif

      jdlimits = [min(jdtt+tbase)-1, max(jdtt+tbase)+1]

      jplephread, efile, info, raw, jdlimits, status=st, errmsg=ee
      if st EQ 0 then message, ee
      jplephinterp, info, raw, jdtt, dpsi, deps, $
        object='NUTATIONS', tbase=tbase

      goto, PRECESS_ANGLES
  endif

  ;; Do this equation in chunks, in case of a large input time array
  dpsi = t*0 & deps = dpsi

  nt = n_elements(t)
  ns = 1000L
  if NOT keyword_set(nonut1) then for i = 0L, nt-1, ns do begin

      ;; Compute indices of input and output arrays
      imax = (i+ns-1)<(nt-1)
      ti = t(i:imax)
      if n_elements(one) NE n_elements(ti) then $
        one = ti*0 + 1

      ;; Compute the fundamental arguments: mean anom of Moon; mean anom
      ;; of Sun; long. of Moon minus long. of Moon's node; mean elongation
      ;; betw. Moon & Sun; mean long. of asc. node of Moon
      ;; ESAA Table 3.222.2
      fundargs = [[poly(ti, farg(*,0))], [poly(ti, farg(*,1))], $
                  [poly(ti, farg(*,2))], [poly(ti, farg(*,3))], $
                  [poly(ti, farg(*,4))]] MOD TWOPI

      ;; ESAA Eqn 3.222-6 (lower equation)
      arg = arg80 # transpose(temporary(fundargs))

      ;; IAU 1980 Nutation in longitude and obliquity (radians)
      ;; ESAA Eqn 3.222-6 (upper equations) and Table 3.222.1
      dpsi(i:imax) = total( (psi80(*,0)#one + psi80(*,1)#ti) * sin(arg), 1)
      deps(i:imax) = total( (eps80(*,0)#one + eps80(*,1)#ti) * cos(arg), 1)
      arg = 0
  endfor

  ;; Above quantities are in mas, convert to radians
  dpsi = dpsi * MAS2R
  deps = deps * MAS2R

  PRECESS_ANGLES:
  ;; Precession from epoch of date to "fixed" epoch
  ;; ESAA Eqn 3.211-2
  t0 = ((fbase-2451545d0) + fepoch)       / 36525d0
  td = ((tbase-fbase) + (jdtt(*)-fepoch)) / 36525d0

  ;; Arguments of ZETA, Z and THETA (part of ESAA Table 3.211.1)
  ;; Symbology note:  ESAA's T is my t0;    ESAA's t is my td
  w1 = poly(t0, [2306.2181d,  1.39656d, -0.000139d])
  w2 = poly(t0, [2004.3109d, -0.85330d, -0.000217d])
  
  ;; IAU 1976 Precession quantities (arcsec)
  ;; Remainder of ESAA Table 3.211.1
  zeta  = (w1 + (( 0.30188D0 - 0.000344D0 * t0 ) + 0.017998D0 * td ) * td )*td
  z     = (w1 + (( 1.09468D0 + 0.000066D0 * t0 ) + 0.018203D0 * td ) * td )*td
  theta = (w2 + ((-0.42665D0 - 0.000217D0 * t0 ) - 0.041833D0 * td ) * td )*td
  ;; ABOVE QUANTITIES ARE IN ARCSEC!
  w1 = 0 & w2 = 0 & td = 0 & t0 = 0  ;; Memory

  ;; Get earth orientation parameters, UT1-UTC
  ;; Convert TT to UT1
  jdtai = jdtt(*) - 32.184d/86400d
  jdutc = jdtai + tai_utc(jdtai + tbase, /invert)/86400d

  ;; Query the EOP database, or set the values to zero
  if keyword_set(useeop) then begin
      eopdata, jdutc, pmx, pmy, ut1_utc, dpsi1, deps1, tbase=tbase

      ;; Adjust the values of the nutation in longitude and obliquity
      deps = deps + deps1
      dpsi = dpsi + dpsi1

      deps1 = 0 & dpsi1 = 0 ;; Memory
  endif else begin
      pmx = 0 & pmy = 0 & ut1_utc = 0 
  endelse
  if keyword_set(no_ut1) then ut1_utc = 0

  ;; Mean obliquity of ecliptic at epoch of date (arcsec)
  ;; ESAA Eqn 3.222-1
  eps0 = poly(t, [84381.448D0, -46.8150D0, -0.00059D0, +0.001813D0]) * AS2R

  ;; True obliquity of ecliptic at epoch of date
  ;; ESAA Eqn 3.222-2
  eps = eps0 + deps

  ;; Greenwich Mean Sidereal Time
  jdut1 = jdutc + ut1_utc/86400d
  jdutc = 0 & jdtai = 0 & & ut1_utc = 0 ;; Memory
  t1 = ((tbase-2451545D) + jdut1) / 36525d

  ;; ESAA Eqn 2.24-1
  gmst1 = poly(t1, [24110.54841d, 8640184.812866d, 0.093104d, -6.2D-6])/86400d
  gmst = (gmst1 MOD 1d) + (jdut1 MOD 1d) + (tbase MOD 1d) + 3.5d
  gmst = TWOPI*( gmst MOD 1d )
  t1 = 0 & gmst1 = 0 ;; Memory

  ;; Equation of the equinoxes - includes two terms depending on the
  ;; mean longitude of the ascending node of the moon.
  ;; ESAA Sec. 3.223 and (IERS Conventions 1996)
  om = poly(t, farg(*,4)) MOD TWOPI
  eqeq = dpsi * cos(eps) + AS2R * (2.64d-3*sin(om) + 6.3d-5*sin(2d*om))
  om = 0 ;; Memory
  
  ;; Greenwich Apparent Sidereal Time
  ;; ESAA Sec. 3.223
  gast = gmst + eqeq

  ;; Units conversions
  case angunits of 
      'ARCSEC': begin
          R2AS = 1d/AS2R ;; Radian to arcsec

          dpsi = dpsi * R2AS
          deps = deps * R2AS
          eps  = eps  * R2AS
          eps0 = eps0 * R2AS
          eqeq = eqeq * R2AS
          pmx  = pmx  * R2AS
          pmy  = pmy  * R2AS
      end
      
      'RADIAN': begin   ;; Arcsec to radians
          zeta = zeta * AS2R 
          z    = z    * AS2R
          theta= theta* AS2R
      end

      ;; Also convert to degrees?

      else: begin
          message, 'ERROR: angular unit '+angunits+$
            ' was not recognized'
      end
  end

  return
end
