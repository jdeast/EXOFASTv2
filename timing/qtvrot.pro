;+
; NAME:
;   QTVROT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Apply quaternion rotation to a 3-vector
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   VNEW = QTVROT(V, Q, [/INVERT])
;
; DESCRIPTION:
;
;   The function QTVROT applies a quaternion rotation (or its inverse)
;   to a 3-vector V to produce a new vector VNEW.
;
;   If both V and VNEW are vector components measured in the same
;   inertial coordinate system, then VNEW returns the components of
;   the vector V rotated by quaternion Q.  I.e., the AXES stay fixed
;   and the VECTOR rotates.  Replace Q by QTINV(Q) in the case of
;   /INVERT.
;
;   If V are components of a vector measured in the "body" coordinate
;   frame, and Q represents the orientation of the body frame
;   w.r.t. the inertial frame, then VNEW are the components of the
;   same vector in the inertial frame.  I.e., the VECTOR stays fixed
;   and the AXES rotate.  For /INVERT, the coordinate transformation
;   is from inertial frame to body frame.
;
;   If either Q is a single quaternion, or V is a single 3-vector,
;   then QTVROT will expand the single to the number of elements of
;   the other operand.  Otherwise, the number of quaternions and
;   vectors must be equal.
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
;
; INPUTS:
;
;  V - array of one or more 3-vectors.  For a single vector, V should
;      be a 3-vector.  For N vectors, V should be a 3xN array.
;
;  Q - array of one or more unit quaternions.  For a single
;      quaternion, Q should be a 4-vector.  For N quaternions, Q
;      should be a 4xN array.
;
;
; RETURNS:
;
;   The resulting rotated vectors.  For single inputs, returns a
;   3-vector.  For N inputs, returns N vectors as a 3xN array.
;
;
; KEYWORD PARAMETERS:
;
;   INVERT - if set, then the antirotation represented by QTINV(Q) is
;            performed.
;
;
; EXAMPLE:
;
;   Q1 = qtcompose([0,0,1],  32d*!dpi/180d)
;   Q2 = qtcompose([1,0,0], 116d*!dpi/180d)
;   Q = qtmult(Q1, Q2)
;
;   V = [[1d,0,0],[0,1,0],[0,0,1]]
;
;   IDL> print, qtvrot(v, q)
;         0.84804810      0.52991926       0.0000000
;         0.23230132     -0.37175982      0.89879405
;         0.47628828     -0.76222058     -0.43837115
;
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Small changes, 28 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;
;  $Id: qtvrot.pro,v 1.7 2002/05/09 23:03:27 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

;; QVROT
;;
;; The FORWARD (default) transform: 
;;
;; * takes a vector vin (components given in inertial coordinates) and
;;  returns the components of the rotated vector vout (components
;;  given in inertial coordinates) -- ie, the AXES stay fixed and the
;;  VECTOR rotates; OR, equivalently,
;;
;; * takes a fixed vector vin (components given in body coordinates)
;;   and returns the components of the vector in inertial coordinates,
;;   where the body system is described by quaternion q -- ie, the
;;   VECTOR stays fixed and the AXES rotate.
;;
;;
;; The INVERSE transform (gotten by setting /INVERT):
;;
;; * takes a vector vin (components given in inertial coordinates) and
;;   returns the components of the anti-rotated vector vout
;;   (components given in inertial coordinates) -- ie, the AXES stay
;;   fixed and the VECTOR rotates.  Anti-rotated here means rotated in
;;   the opposite direction of q; OR, equivalently,
;;
;; * takes a fixed vector vin (components given in inertial
;;   coordinates) and returns the components of the vector in body
;;   coordinates, where the body system is described by quaternion q
;;   -- ie, the VECTOR stays fixed and the AXES rotate.
;;

function qtvrot, vin, q, invert=invert

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'VNEW = QTVROT(V, Q)', info=info
      return, 0
  endif
  nq = n_elements(q)/4
  nv = n_elements(vin)/3
  if nq LT 1 OR nv LT 1 then goto, USAGE

  if n_elements(q) GT 4 AND n_elements(vin) GT 3 then begin
      if n_elements(q)/4 NE n_elements(vin)/3 then begin
          message, 'ERROR: incompatible number of quaternions & vectors'
          return, -1L
      end
      vout = vin*q(0)*0.
      nq = n_elements(q)/4
      nv = nq
  endif else if n_elements(q) GT 4 then begin
      nq = n_elements(q)/4
      nv = 1L
      vout = vin(*) # (fltarr(nq)+1) * q(0)*0.
  endif else begin
      nq = 1L
      nv = n_elements(vin)/3
      vout = vin*q(0)*0.
  endelse
  vout = reform(vout, 3, max([nv,nq]), /overwrite)
  q1 = q(0,*) & q2 = q(1,*) & q3 = q(2,*) & q4 = q(3,*)
  if n_elements(q1) EQ 1 then begin
      q1 = q1(0)  & q2 = q2(0)  & q3 = q3(0)  & q4 = q4(0)
  endif else begin
      q1 = q1(*)  & q2 = q2(*)  & q3 = q3(*)  & q4 = q4(*)
  endelse
  v0 = vin(0,*) & v1 = vin(1,*) & v2 = vin(2,*)
  if n_elements(v0) EQ 1 then begin
      v0 = v0(0)  & v1 = v1(0)  & v2 = v2(0)
  endif else begin
      v0 = v0(*)  & v1 = v1(*)  & v2 = v2(*)
  endelse

  if NOT keyword_set(INVERT) then begin

      ;; FORWARD TRANSFORMATION
      VOUT(0,*)=((Q1*Q1-Q2*Q2-Q3*Q3+Q4*Q4)*V0 $
                 + 2.D0*(Q1*Q2-Q3*Q4)*V1 $
                 + 2.D0*(Q1*Q3+Q2*Q4)*V2)
      VOUT(1,*)=(2.D0*(Q1*Q2+Q3*Q4)*V0 $
                 + (-Q1*Q1+Q2*Q2-Q3*Q3+Q4*Q4)*V1 $
                 + 2.D0*(Q2*Q3-Q1*Q4)*V2)
      VOUT(2,*)=(2.D0*(Q1*Q3-Q2*Q4)*V0 $
                 + 2.D0*(Q2*Q3+Q1*Q4)*V1 $
                 + (-Q1*Q1-Q2*Q2+Q3*Q3+Q4*Q4)*V2)
  endif else begin
      ;; INVERSE TRANSFORMATION
      VOUT(0,*)=((Q1*Q1-Q2*Q2-Q3*Q3+Q4*Q4)*V0 $
                 + 2.D0*(Q1*Q2+Q3*Q4)*V1 $
                 + 2.D0*(Q1*Q3-Q2*Q4)*V2)
      VOUT(1,*)=(2.D0*(Q1*Q2-Q3*Q4)*V0 $
               + (-Q1*Q1+Q2*Q2-Q3*Q3+Q4*Q4)*V1 $
               + 2.D0*(Q2*Q3+Q1*Q4)*V2)
      VOUT(2,*)=(2.D0*(Q1*Q3+Q2*Q4)*V0 $
               + 2.D0*(Q2*Q3-Q1*Q4)*V1 $
               + (-Q1*Q1-Q2*Q2+Q3*Q3+Q4*Q4)*V2)
  endelse
  vout = vout

  return, vout
end
