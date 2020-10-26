; PURPOSE:
;   This generates the lines for the EXOFASTv2 prior file for wide
;   uniform LD priors based on an extrapolation of the Claret LD
;   tables. This is designed for low mass stars where the Claret
;   tables are not applicable, but should be in the right ballpark.


pro ldpriors, teff, feh, bands,logg=logg, mstar=mstar, rstar=rstar

if n_elements(logg) eq 0 then begin
   if n_elements(mstar) eq 0 or n_elements(rstar) eq 0 then begin
      message, 'logg or mstar and rstar must be specified'
   endif
   constants = mkconstants()
   logg = alog10(mstar/(rstar^2)*constants.gravitysun) ;; cgs
endif

nbands = n_elements(bands)

for i=0L, nbands-1 do begin
   
   coeffs = quadld(logg, teff, feh, bands[i],/skipbound)

   print, '# ' + bands[i] 
   print, 'u1_' + strtrim(i,2) + ' ' + strtrim(coeffs[0],2) + ' -1 ' + strtrim(coeffs[0]-0.1d0,2) + ' ' +  strtrim(coeffs[0]+0.1d0,2) 
   print, 'u2_' + strtrim(i,2) + ' ' + strtrim(coeffs[1],2) + ' -1 ' + strtrim(coeffs[1]-0.1d0,2) + ' ' +  strtrim(coeffs[1]+0.1d0,2) 

endfor

end
