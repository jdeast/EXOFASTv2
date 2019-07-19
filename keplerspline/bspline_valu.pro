;+
; NAME:
;   bspline_valu
;
; PURPOSE:
;      1) Evaluate a bspline set (see create_bsplineset) at specified
;            x and x2 arrays
;
; CALLING SEQUENCE:
;   
;   yfit  = bspline_valu( x, sset, x2=x2, action=action, upper=upper, 
;               lower=lower, mask=mask)
;
; INPUTS:
;   x          - independent variable
;   sset       - Structure to be returned with all fit parameters
;
; RETURNS:
;   yfit       - Evaluated b-spline fit
;
; OPTIONAL KEYWORDS:
;   x2         - Orthogonal dependent variable for 2d fits
;   action     - This keyword is overwritten with b-spline action matrix
;   lower,upper- Internal keywords used by action, maybe should replace
;                action with a structure including lower and upper
;
; OPTIONAL OUTPUTS:
;   mask       - a mask array of good (1's) bspline evalutions
; 
; COMMENTS:
;   the mask attempts to show regions where the bspline was ill-defined
;    and breakpoints had been dropped.
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   11-Sep-2000 Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function bspline_valu, x, sset, x2=x2, action=action, upper=upper, $
    lower=lower, mask=mask

      nx = n_elements(x)
      mask = lonarr(nx) 

      if size(sset,/tname) NE 'STRUCT' then begin
         message, 'Please send in a proper B-spline structure', /continue
         return, x*0.0
      endif

      xsort = sort(x)
      npoly = 1L
      xwork = x[xsort]

      if keyword_set(x2) then begin
        if ((where(tag_names(sset) EQ 'NPOLY'))[0] NE -1) then $
           npoly=sset.npoly
        x2work = x2[xsort]
      endif else x2work = 0

      if NOT keyword_set(action) then $
           action = bspline_action(xwork, sset, x2=x2work, upper=upper, $
                       lower=lower)
 
      yfit = x * 0.0
      nord = sset.nord
      bw = npoly * nord


      spot = lindgen(bw)
      goodbk = where(sset.bkmask NE 0, nbkpt) 
      coeffbk = where(sset.bkmask[nord:*] NE 0) 
      n = nbkpt - nord

      
      sc = size(sset.coeff)
      if sc[0] EQ 2 then goodcoeff = sset.coeff[*,coeffbk] $
      else goodcoeff = sset.coeff[coeffbk]

      maskthis = xwork * 0.0 
   

      for i= 0L, n - nord do begin

         ict = upper[i] - lower[i] + 1

          if (ict GT 0) then begin
             yfit[lower[i]:upper[i]] = action[lower[i]:upper[i],*] # $
                  goodcoeff[i*npoly+spot]
          endif

      endfor

      yy = yfit
      yy[xsort] = yfit 

      mask[*] = 1
      gb = sset.fullbkpt[goodbk]

      outside = where(x LT gb[nord-1] OR x GT gb[n])
      if outside[0] NE -1 then mask[outside] = 0
   
      hmm = where(goodbk[1:*] - goodbk GT 2, nhmm) 
      for jj=0, nhmm - 1 do begin
        inside = where(x GE sset.fullbkpt[goodbk[hmm[jj]]] $
                  AND  x LE sset.fullbkpt[goodbk[hmm[jj]+1]-1])
        if inside[0] NE -1 then mask[inside] = 0
      endfor
        
      return, yy
end

