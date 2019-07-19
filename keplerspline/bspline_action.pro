;+
; NAME:
;   bspline_action
;
; PURPOSE:
;      1) Construct banded bspline matrix, with dimensions [ndata, bandwidth]
;
; CALLING SEQUENCE:
;   
;   action = bspline_action( x, sset, x2=x2, lower=lower, upper=upper)
;
; INPUTS:
;   x          - independent variable
;   sset       - Structure to be returned with all fit parameters
;
; RETURNS:
;   action     - b-spline action matrix
;
; OPTIONAL KEYWORDS:
;   x2         - Orthogonal dependent variable for 2d fits
;
; OPTIONAL OUTPUTS:
;   lower      - A list of pixel positions, each corresponding to the first
;                occurence of position greater than breakpoint indx
;   upper      - Same as lower, but denotes the upper pixel positions
;
; COMMENTS:
;   Does not yet support the slatec function to directly return
;   derivatives of the b-spline (ideriv). 
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Sep-2000 Written by Scott Burles, FNAL
;    3-Jul-2001 Fundamental array organization bug fixed, S. Burles 
;-
;------------------------------------------------------------------------------
function bspline_action, x, sset, x2=x2, lower=lower, upper=upper

      if size(sset,/tname) NE 'STRUCT' then begin
         message, 'Please send in a proper B-spline structure', /continue
         return, -1L
      endif
 
      npoly = 1L
      nx = n_elements(x)

;
;	Check for the existence of x2 
;
      if keyword_set(x2) then begin
        if n_elements(x2) NE nx then begin
          print, 'dimensions do not match between x and x2'
          return, -1
        endif

        if ((where(tag_names(sset) EQ 'NPOLY'))[0] NE -1) then npoly=sset.npoly
      endif
 

      nord = sset.nord
      goodbk = where(sset.bkmask NE 0, nbkpt)
      if nbkpt LT 2*nord then return, -2L
      n = nbkpt - nord

      gb = sset.fullbkpt[goodbk]

      bw = npoly * nord   
      action = x # replicate(0,bw) 

      lower = lonarr(n-nord+1)
      upper = lonarr(n-nord+1) - 1

      indx = intrv(x, gb, nord)

      bf1 = bsplvn(gb, nord, x, indx)
      action = bf1

      ;--------------------------------------------------------------
      ;  sneaky way to calculate upper and lower indices when 
      ;   x is sorted
      ;
      aa = uniq(indx)
      upper[indx[aa]-nord+1] = aa

      rindx = reverse(indx)
      bb = uniq(rindx)
      lower[rindx[bb]-nord+1] = nx - bb - 1

      ;---------------------------------------------------------------
      ;  just attempt this if 2d fit is required
      ;
      if keyword_set(x2) then begin


         x2norm = 2.0 * (x2[*] - sset.xmin) / (sset.xmax - sset.xmin) - 1.0
         CASE sset.funcname OF
           'poly' : begin
                   temppoly = (x2norm*0.0 + 1.0) # replicate(1,npoly)
                   for i=1,npoly-1 do temppoly[*,i] = temppoly[*,i-1] * x2norm
                  end
           'poly1' : begin
                   temppoly = x2norm # replicate(1,npoly)  
                   for i=1,npoly-1 do temppoly[*,i] = temppoly[*,i-1] * x2norm
                  end
           'chebyshev' : temppoly = fchebyshev(x2norm, npoly)
           'legendre'  : temppoly = flegendre(x2norm, npoly)
           else :        temppoly = flegendre(x2norm, npoly)
         ENDCASE

         action = fltarr(nx,bw)
         counter= -1L
         for ii=0,nord-1 do begin
           for jj=0,npoly-1 do begin
             counter = counter +1
             action[*, counter] = bf1[*,ii] * temppoly[*,jj]
           endfor
         endfor
      endif

      return, action
end

