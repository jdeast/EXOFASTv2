;+
; NAME:
;   QTEULER
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Compose a series of euler-type rotations into a single quaternion
;
; MAJOR TOPICS:
;   Geometry
;
; CALLING SEQUENCE:
;   Q = QTEULER(AXES, ANG0, ANG1, ... )
;
; DESCRIPTION:
;
;  The function QTEULER composes a series of Euler-type rotations into
;  a single set of quaternion representations.
;
;  The user specifies a set of axes, and the angles to rotation about
;  those axes, and QTEULER constructs the corresponding quaternion.
;
;  There must be a one-to-one correspondence between the elements of
;  AXES and the number of rotations.  AXES specifies the rotation axes
;  as an string, which must be one of 'X', 'Y', or 'Z'.  Other axes
;  are invalid.  For example, the following call:
;
;    QTEULER(['X','Z'], THETA, PHI)
;
;  will rotate first about the *Z* axis by the angle PHI, and then
;  around the *resulting X* axis by angle THETA.
;
;  Several things are worth noting here.  First, rotations are applied
;  first from the right, not the left.  This conforms to the usual
;  matrix notation for applying rotations to a vector on the right
;  hand side.  For example, in matrix notation,
;
;      XNEW = A3 A2 A1 XOLD
;
;  applies first A1, then A2 and finally A3 to the XOLD vector,
;  resulting in the new vector XNEW.  The same semantics apply here.
;
;  A second thing to bear in mind is that the axes themselves change
;  during the rotations.  Thus, the coordinates specified in AXES
;  should be considered attached to the "body" and not the inertial
;  frame.
;
;
; INPUTS:
;
;  AXES - a string array, specifies the rotation axes.  Rotations are
;         applied last element first.  Each element of AXES must be
;         one of 'X', 'Y' or 'Z'.
;
;  ANG0, ..., ANGi - the successive rotation angles.  Angle ANGi
;         corresponds to axis AXES(i).
;
;         If ANGi is a scalar, then it will be promoted to a vector
;         the same size as the other rotation angles being performed.
;         Otherwise, if the angles ANGi are vectors, then they must
;         all be of the same size.
;
; RETURNS:
;
;  The resulting quaternion (or, if ANGi are vectors, array of
;  quaternions), which represent the requested rotations.
;
;
; KEYWORD PARAMETERS:
;
;  NONE
;
; EXAMPLE:
;
;    ;;              Precession        Nutation    
;    qtot = qteuler(['z','y','z',      'x','z','x'        ], $
;                    -zeta, +theta, -z, +eps0, -dpsi, -eps)
;
;   Applies a series of rotations to correct for earth nutation and
;   precession.  The order of rotations on a vector would be
;   X-Z-X-Z-Y-Z (i.e., the reverse order printed).
;
; SEE ALSO
;   QTANG, QTAXIS, QTCOMPOSE, QTERP, QTEXP, QTFIND, QTINV, QTLOG,
;   QTMAT, QTMULT, QTPOW, QTVROT
;
; MODIFICATION HISTORY:
;   Written, 27 Jan 2002, CM
;   More error checking, 03 Mar 2002, CM
;
;  $Id: qteuler.pro,v 1.4 2002/05/09 23:03:27 craigm Exp $
;
;-
; Copyright (C) 2002, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

;; Extract axis ei and angle angi
pro qteuler_extract, ax, i, ei, angi, $
                     ang0, ang1, ang2, ang3, ang4, $
                     ang5, ang6, ang7, ang8, ang9, $
                     status=status, errmsg=errmsg
                     
  status = 0
  zero = ang0(0)*0
  ex = [1,zero,zero] & ey = [zero,1,zero] & ez = [zero,zero,1]

  ei = [0D, 0D, 0D]

  if execute('ei = e'+ax(i)+' & angi = ang'+strtrim(i,2)) NE 1 then begin
      stop
      errmsg = 'Invalid axis specification'
      return
  endif

  status = 1
  return
end

function qteuler, axes, block=block, $
            ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9, $
            ang10, ang11, ang12, ang13, ang14, ang15

  if n_params() EQ 0 then begin
      info = 1
      USAGE_ERR:
      message, 'USAGE: Q = QTEULER(AXES, ANG0, ...)', /info
      message, '  AXES = ["X",...]    ("X" or "Y" or "Z")', /info
      message, '  ANGn = rotation angle (radians)', info=info
      return, 0
  endif
  if n_elements(axes) LT 1 OR n_elements(ang0) LT 1 then $
    goto, USAGE_ERR
  nang = n_params()-1

  ;; Check to be sure each axis label is 'X' 'Y' or 'Z'
  ax = strupcase(strmid(strtrim(axes,2),0,1))
  wh = where(ax NE 'X' AND ax NE 'Y' AND ax NE 'Z', ct)
  if ct GT 0 then begin
      errmsg = 'AXES must be one of "X", "Y" or "Z"'
      goto, BAD_AXIS
  endif
  if n_elements(ax) NE nang then begin
      errmsg = 'Number of AXES and rotations ANGi must agree'
      goto, BAD_AXIS
  endif

  qteuler_extract, ax, 0, ev, angv, status=status, errmsg=errmsg, $
    ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9

  if status EQ 0 then begin
      BAD_AXIS:
      message, 'ERROR: '+errmsg, /info
      goto, USAGE_ERR
  endif

  qq = qtcompose(ev, angv)
  for i = 1, nang-1 do begin
      qteuler_extract, ax, i, ev, angv, status=status, errmsg=errmsg, $
        ang0, ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9
      if status EQ 0 then goto, BAD_AXIS
      
      qq = qtmult(qq, qtcompose(ev, angv))
  endfor

  return, qq
end
