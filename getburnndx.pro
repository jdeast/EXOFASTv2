;+
; NAME:
;   GETBURNNDX
; PURPOSE:
;   Returns the index of the first useable link in a Markov chain.
;
; DESCRIPTION: 
;   Determines the point at which all chains have had at least one
;   chi^2 lower than the median chi^2 of the best chain. If some
;   chains never reach this point, they are ignored in the
;   calculation. 
;   
;   nsteps*0.1 <= burnndx <= nsteps*0.9
; 
;   That is...
;
;   It will always discard at least 10% of the steps in case the
;   initialization of the chains is insufficiently scattered around
;   the best fit, which would bias the results.
; 
;   It will never discard more than 90% of the steps to enable the
;   robust calculation of mixing criteria.
;
; CALLING SEQUENCE: 
;   burnndx = getburnndx(chi2 [,GOODCHAINS=,BADCHAINS=])
;
; INPUTS:
;   
;   CHI2 - An NSTEPS x NCHAINS array of chi2 values
; 
; OPTIONAL OUTPUTS:
; 
;   GOODCHAINS - An array of indices of chains that are good, meaning
;                they have at least one point below the median of the
;                best chain.  
;   BADCHAINS  - An array of indices of chains that are bad, meaning
;                their chi2 is never below the median of the best
;                chain
;
;   MEDCHI2    - The median chi2 of the best chain
;
; REVISION HISTORY:
;   2018/03 - Written -- Jason Eastman (CfA)
;-
function getburnndx, chi2, goodchains=goodchains, badchains=badchains, medchi2=medchi2

sz = size(chi2)
nsteps = sz[1]
nchains = sz[2]
;; error checking
if sz[0] ne 2 then message, 'CHI2 must be an NSTEPS x NCHAINS (2D) array'
if nsteps lt 5 then message, 'CHI2 must have at least 5 steps'
if nchains lt 3 then message, 'CHI2 must have at least 3 chains'

burnndx = round(0.1*nsteps) ;; discard at least the first 10% of the chains

;; the minimum chi2 of each chain
minchi2 = min(chi2[burnndx:nsteps-1,*],dimension=1)

;; some chains get stuck (why??) -- identify them and discard them
uniqvals = dblarr(nchains)
for i=0, nchains-1 do uniqvals[i]= n_elements(uniq(chi2[*,i]))
djs_iterstat,uniqvals,sigrej=5,mask=mask
goodchains = where(mask eq 1)
if goodchains[0] eq -1 then goodchains = lindgen(nchains)

;; the median chi2 of the best chain
medchi2 = min(median(chi2[burnndx:nsteps-1,goodchains],dimension=1))

;; good chains must have some values below this median
goodchains = where(minchi2 lt medchi2 and mask eq 1, ngood, complement=badchains)

;; if less than 3 chains are good, use them all 
;; this isn't going to go well...
if ngood lt 3 then begin
   goodchains = lindgen(nchains)
   ngood = nchains
endif

;; the burnndx is the first link where all good chains have been below
;; the median chi2 at least once
for j=0L, ngood-1 do begin
   tmpndx = (where(chi2[*,goodchains[j]] lt medchi2))(0)
   if tmpndx gt burnndx then burnndx = tmpndx
endfor

;; don't throw away more than 90% of the chain
;; this also isn't going to go well...
burnndx = burnndx < round(nsteps*0.9d0)

return, burnndx

end
