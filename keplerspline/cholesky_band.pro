
function cholesky_band, lower, mininf=mininf, verbose=verbose


    if NOT keyword_set(mininf) then mininf=0.0
    ; compute cholesky decomposition of banded matrix
    ;   lower[bandwidth, n]  n is the number of linear equations
 
    ; I'm doing lower cholesky decomposition from lapack, spbtf2.f

    bw = (size(lower))[1]
    n = (size(lower))[2] - bw

    kd = bw - 1


    negative = where(lower[0,0:n-1] LE mininf)
    if negative[0] NE -1 then begin
       if keyword_set(verbose) then begin
           message,'bad entries',/continue
           print, negative
       endif

       return, negative
    endif

    kn = bw - 1 
    spot = 1 + lindgen(kn)
    bi = lindgen(kn)
    for i=1L, kn-1 do bi = [bi,lindgen(kn-i)+(kn+1)*i]


    for j=0L, n-1 do begin
         lower[0,j] = sqrt(lower[0,j])
         lower[spot,j] = lower[spot,j] / lower[0,j]
         x = lower[spot,j]

         if (where(finite(x) EQ 0))[0] NE -1 then begin
 	    if keyword_set(verbose) then $
               message, 'NaN found in cholesky_band', /continue
            return, j
         endif

         hmm = x # transpose(x)
         here = bi+(j+1)*bw
         lower[here] = lower[here] - hmm[bi]
    endfor

  return,-1L
end
