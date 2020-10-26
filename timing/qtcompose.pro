;+
; NAME:
;   QTCOMPOSE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Convert a rotation angle and axis into quaternion
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   Q = QTCOMPOSE(VAXIS, PHI)
;
; DESCRIPTION:
;
;  The function QTCOMPOSE accepts a unit vector rotation axis VAXIS
;  and a rotation angle PHI, and returns the corresponding quaternion.
;
;  The user must take care to pass the same number of axes as rotation
;  angles.
;
;  Use QTAXIS and QTANG to extract the properties of an existing
;  quaternion.  Use QTCOMPOSE to combine a rotation axis and angle
;  into a new quaternion.
;
;  Conventions for storing quaternions vary in the literature and from
;  library to library.  This library uses the convention that the
;  first three components of each quaternion are the 3-vector axis of
;  rotation, and the 4th component is the rotation angle.  Expressed
;  in formulae, a single quaternion is given by:
;
;     Q(0:2) = [VX, VY, VZ]*SIN(PHI/2)
;     Q(3)   =              COS(PHI/2)
;
;  where PHI is the rotation angle, and VAXIS = [VX, VY, VZ] is the
;  rotation eigen axis expressed as a unit vector.  This library
;  accepts quaternions of both signs, but by preference returns
;  quaternions with a positive 4th component.
;
; INPUTS:
;
;  VAXIS - array of one or more unit vectors specifying the rotation
;          axes.  For a single rotation, VAXIS should be a 3-vector.
;          For N vectors, VAXIS should be a 3xN array.
;
;  PHI - one or more rotation angles, in radians.  For a single
;        rotation, PHI should be a scalar.  For N rotations, PHI
;        should be an N-vector.
;
; RETURNS:
;
;  For a single rotation, returns a quaternion as a 4-vector.  For N
;  rotations, returns a 4xN vector of quaternions.
;
;
; KEYWORD PARAMETERS:
;
;  NONE
;
; EXAMPLE:
;
;   IDL> print, qtcompose([0d,1,0], !dpi/4)
;          0.0000000      0.38268343       0.0000000      0.92387953
;
;   Prints the quaternion composed of a rotation of !dpi/4 radians
;   around the axis [0,1,0]
;
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Allow output to be DOUBLE, 27 Jan 2002, CM
;   Allow vector vs scalar arguments, 28 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;
;  $Id: qtcompose.pro,v 1.11 2008/12/14 20:00:31 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
function qtcompose, axis, phi

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'Q = QTCOMPOSE(AXIS, PHI)', info=1
      return, 0
  endif

  nph = n_elements(phi)
  nv  = n_elements(axis)/3
  if nph LT 1 OR nv LT 1 then goto, USAGE

  nq  = nv > nph
  q = make_array(value=axis(0)*phi(0)*0., 4,nq)

  sph = sin(phi/2) & cph = cos(phi/2)
  if nph EQ 1 AND nv EQ 1 then return, [ axis(0:2) * sph(0), cph(0) ]
  if nph GT 1 AND nv EQ 1 then begin
      ;; Single axis, multiple rotation angles
      q(0,*) = axis(0)*sph
      q(1,*) = axis(1)*sph
      q(2,*) = axis(2)*sph
      q(3,*) = cph
  endif else if nph EQ 1 AND nv GT 1 then begin
      ;; Multiple axis, single rotation
      q(0:2,*) = axis*sph(0)
      q(3,*)   = cph(0)
  endif else if nph EQ nv then begin
      ;; Multiple axes, multiple rotations
      q(0:2,*) = axis*rebin(reform(temporary(sph),1,nq),3,nq)
      q(3,*)   = temporary(cph)
  endif else begin
      message, 'ERROR: number of axes and angles do not match'
  endelse

  return, q
end
