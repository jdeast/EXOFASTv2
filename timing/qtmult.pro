;+
; NAME:
;   QTMULT
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Multiply quaternions
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   QRESULT = QTMULT(Q1, [/INV1,]   Q2, [/INV2])
;
; DESCRIPTION:
;
;   The function QTMULT performs multiplication of quaternions.
;   Quaternion multiplication is not component-by-component, but
;   rather represents the composition of two rotations, namely Q2
;   followed by Q1.
;
;   More than one multiplication can be performed at one time if Q1
;   and Q2 are 4xN arrays.  In that case both input arrays must be of
;   the same dimension.
;
;   If INV1 is set, then the inverse of Q1 is used.  This is a
;   convenience, to avoid the call QTINV(Q1).  Of course, INV2 can
;   be set to use the inverse of Q2.
;
;   Note that quaternion multiplication is not commutative.
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
;  Q1 - array of one or more unit quaternions, the first operand in
;       the multiplication.  For a single quaternion, Q1 should be a
;       4-vector.  For N quaternions, Q1 should be a 4xN array.
;       If INV1 is set, then the inverse of Q1 is used.
;
;  Q2 - same as Q1, for the second operand.
;       If INV2 is set, then the inverse of Q2 is used.
;
; RETURNS:
;
;   The resulting multiplied unit quaternions.  For a single inputs,
;   returns a 4-vector.  For N input quaternions, returns N
;   quaternions as a 4xN array.
;
;
; KEYWORD PARAMETERS:
;
;  INV1 - if set, use QTINV(Q1) in place of Q1.
;
;  INV2 - if set, use QTINV(Q2) in place of Q2.
;
; EXAMPLE:
;
;   Q1 = qtcompose([0,0,1],  32d*!dpi/180d)
;   Q2 = qtcompose([1,0,0], 116d*!dpi/180d)
;
;   IDL> print, qtmult(q1, q2)
;        0.81519615      0.23375373      0.14606554      0.50939109
;
;   Form a rotation quaternion of 32 degrees around the Z axis, and 
;   116 degrees around the X axis, then multiply the two quaternions.
;   
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTMULTN, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, July 2001, CM
;   Documented, Dec 2001, CM
;   Documentation, allow 1xN or Nx1 multiplies, 27 Jan 2002, CM
;   Usage message, error checking, 15 Mar 2002, CM
;   Add the INV1 and INV2 keywords, 30 Aug 2007, CM
;
;  $Id: qtmult.pro,v 1.8 2007/09/03 07:18:25 craigm Exp $
;
;-
; Copyright (C) 2001, 2002, 2007, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

function qtmult, aqt, bqt, inv1=inverse1, inv2=inverse2

; THIS ROUTINE MULTIPLIES QUATERNIONS
; CQT CORRESPONDS TO THE ROTATION AQT FOLLOWED BY BQT
; ASSUMING S/C COORDINATES ARE INITIALLY ALIGN WITH INERTIAL COORD.
; THEN ROTATION AQT DESCRIBES ROTATION SUCH THAT THE SUBROUTINE
;   QTXRA GIVES THE INERTIAL COORDINATES OF THE S/C X-AXIS
;   THE FIRST 3 COMPONENTS OF AQT GIVE THE EIGENAXIS EXPRESSED
;   IN S/C COORDINATES BEFORE THE ROTATION (=INTERTIAL COORD.).
; THE BQT ROTATION FOLLOWS THE AQT ROTATION. CQT THEN DESCRIBES
;   THIS COMBINATION SUCH THAT QTXRA GIVES THE INERTIAL COORDINATES
;   OF THE S/C X-AXIS AFTER BOTH ROTATIONS. 
;   THE FIRST 3 COMPONENTS OF BQT GIVE THE EIGENAXIS EXPRESSED
;   IN S/C COORDINATES AFTER THE AQT ROTATION.

  if n_params() EQ 0 then begin
      info = 1
      USAGE:
      message, 'USAGE:', /info
      message, 'QNEW = QTMULT(Q1, Q2)', info=info
      return, 0
  endif

  sz1 = size(aqt)
  sz2 = size(bqt)
  if sz1(0) LT 1 OR sz2(0) LT 1 then $
    message, 'ERROR: Q1 and Q2 must be quaternions'
  if sz1(1) NE 4 OR sz2(1) NE 4 then $
    message, 'ERROR: Q1 and Q2 must be quaternions'
  n1 = n_elements(aqt)/4
  n2 = n_elements(bqt)/4
  if n1 NE n2 AND n1 NE 1 AND n2 NE 1 then $
    message, 'ERROR: Q1 and Q2 must both have the same number of quaternions'

  nq = n1>n2
  cqt = make_array(value=aqt(0)*bqt(0)*0, dimension=[4,nq])

  if n1 GT 1 then begin
      aqt0 = aqt(0,*) & aqt1 = aqt(1,*) & aqt2 = aqt(2,*) & aqt3 = aqt(3,*)
  endif else begin
      aqt0 = aqt(0) & aqt1 = aqt(1) & aqt2 = aqt(2) & aqt3 = aqt(3)
  endelse
  if n2 GT 1 then begin
      bqt0 = bqt(0,*) & bqt1 = bqt(1,*) & bqt2 = bqt(2,*) & bqt3 = bqt(3,*)
  endif else begin
      bqt0 = bqt(0) & bqt1 = bqt(1) & bqt2 = bqt(2) & bqt3 = bqt(3)
  endelse
  if keyword_set(inverse1) then begin
      aqt0 = -aqt0 & aqt1 = -aqt1 & aqt2 = -aqt2
  endif
  if keyword_set(inverse2) then begin
      bqt0 = -bqt0 & bqt1 = -bqt1 & bqt2 = -bqt2
  endif

  CQT(0,0) = AQT0*BQT3 + AQT1*BQT2 - AQT2*BQT1 + AQT3*BQT0
  CQT(1,0) =-AQT0*BQT2 + AQT1*BQT3 + AQT2*BQT0 + AQT3*BQT1
  CQT(2,0) = AQT0*BQT1 - AQT1*BQT0 + AQT2*BQT3 + AQT3*BQT2
  CQT(3,0) =-AQT0*BQT0 - AQT1*BQT1 - AQT2*BQT2 + AQT3*BQT3
  
  return, cqt
end
