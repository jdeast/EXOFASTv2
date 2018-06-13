;+
; NAME:
;  angsep
; PURPOSE:
;  Compute the angular distance between two spherical coordinates.
; DESCRIPTION:
;
; CATEGORY:
;  Mathematical
; CALLING SEQUENCE:
;  ans=angsep(ra1,dec1,ra2,dec2)
;
; INPUTS:
;  ra1  - RA of first position (radians)
;  dec1 - Dec of first position (radians)
;  ra2  - RA of second position (radians)
;  dec2 - Dec of second position (radians)
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;
; OUTPUTS:
;  return value is the angular distance in radians.
;
; KEYWORD OUTPUT PARAMETERS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;  Written by Marc W. buie, Lowell Observatory, 1997/09/08
;  2009/02/26, MWB, added protection against round off error generating NaN
;-
function angsep,ra1,dec1,ra2,dec2

   arg = sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2)
   arg = ( (arg<1) > (-1) )

   angsep=acos(arg)

   return,angsep

end
