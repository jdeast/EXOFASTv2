;+
; NAME:
;   bspline_maskpoints
;
; PURPOSE:
;   Perform simple logic of which breakpoints to mask
;
; CALLING SEQUENCE:
;   error_code = bspline_maskpoints(sset, errb, npoly)
;   
; INPUTS:
;   sset      - Bspline structure
;   errb      - Errors returned from cholesky routines
;   npoly     - Polynomial norder of 2d fit (default 1)
;
; RETURNS:
;   error_code =  -1:   Breakpoints have been dropped, try again
;                 -2:   Solution not possible, abort
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Oct-2000 Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function bspline_maskpoints, sset, errb, npoly

   if NOT keyword_set(npoly) then npoly=1L

   goodbk = where(sset.bkmask NE 0, nbkpt)
   nord = sset.nord

   if nbkpt LE 2*nord then return, -2L

   hmm = errb[uniq(errb/npoly)]/npoly  
   n = nbkpt - nord

   if (where(hmm GE n))[0] NE -1 then return, -2L

   test = lonarr(nbkpt) 
   for jj=-ceil(nord/2.0),(nord/2.0)-1 do begin
      inside = (((hmm + jj) > 0) + nord) < (n-1)
      test[inside] = 1
   endfor
 
   maskthese = where(test EQ 1)


   if maskthese[0] EQ -1 then return, -2L

   reality = goodbk[maskthese]
   if total(sset.bkmask[reality]) EQ 0 then return, -2L

   sset.bkmask[reality] = 0 
   return, -1L

end
