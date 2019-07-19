;+
; NAME:
;   intrv
;
; PURPOSE:
;   Find the segment between breakpoints which contain each value in
;   the array x.  The minimum breakpoint is nbkptord -1, and the maximum
;   is nbkpt - nbkptord - 1.  This routine is required by the bspline
;   IDL routines, and is similar in function to the slatec version.
;
; CALLING SEQUENCE:
;   
;   indx  = intrv(x, fullbkpt, nbkptord)
;
; INPUTS:
;   x          - data x values
;   fullbkpt   - Breakpoint vector returned by efc
;   nbkptord   - Order of b-splines (4 is cubic)
;
; RETURNS:
;   indx       - position of array elements with respect to breakpoints.
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   does the same function as intrv, although slower but easier to follow
;    sorting is done here
;   Again, assumes x is monotonically increasing
;
; EXAMPLES:
;
; REVISION HISTORY:
;   31-Aug-2000 Written by Scott Burles, FNAL
;-
;------------------------------------------------------------------------------
function intrv, x, fullbkpt, nbkptord 
    
      nx = n_elements(x)
      nbkpt= n_elements(fullbkpt)
      n = (nbkpt - nbkptord)

      indx = lonarr(nx)

      ileft = nbkptord - 1L
      for i=0L, nx-1 do begin
        while (x[i] GT fullbkpt[ileft+1] AND ileft LT n-1 ) do $
            ileft = ileft + 1L
        indx[i] = ileft
      endfor
      indxold = indx


     ; here's another sneaky attempt, which takes way too long for
     ; superflat...

     ; fullist = [fullbkpt, x]
     ; hmm = sort(fullist)
     ; back = lonarr(n_elements(hmm))
     ; back[hmm] = lindgen(n_elements(hmm))
     ; for i=0,n-1 do indx[back[i]-i:nx-1] = i  

      
     return, indx
end 
