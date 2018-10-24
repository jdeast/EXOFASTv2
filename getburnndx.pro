;+
; NAME:
;   GETBURNNDX
; PURPOSE:
;   Returns the index of the first useable link in a Markov chain.
;
; DESCRIPTION: 
;   Determines the point at which all chains have had at least one
;   chi^2 lower than the median chi^2 of the best chain. If some
;   chains never reach this point or reach this point very late, they
;   are ignored in the calculation.
;   
;   nsteps*0.1 <= burnndx <= nsteps*0.9
; 
;   That is...
;
;   It will always discard at least 10% of the steps in case the
;   initialization of the chains is insufficiently scattered around
;   the best fit, which would bias the results toward the starting
;   value.
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
;   2018/04 - Modified goodchains to maximize the total number of
;             saved links
;-
function getburnndx, chi2, goodchains=goodchains, badchains=badchains, medchi2=medchi2, logname=logname

sz = size(chi2)
nsteps = sz[1]
nchains = sz[2]
;; error checking
if sz[0] ne 2 then message, 'CHI2 must be an NSTEPS x NCHAINS (2D) array'
if nsteps lt 5 then message, 'CHI2 must have at least 5 steps'
if nchains lt 3 then message, 'CHI2 must have at least 3 chains'

;; the preliminary burnndx is the index of the best chi2 among all chains
bestchi2 = !values.d_infinity
for i=0L, nchains-1 do begin
   tmpchi2 = min(chi2[*,i], ndx)
   if tmpchi2 lt bestchi2 then begin
      bestchi2 = tmpchi2
      bestchain = i
      burnndx = ndx
   endif
endfor
burnndx = (burnndx > round(0.1*nsteps)) < round(0.75*nsteps)

;; now find the median chi2 of that chain after the minimum 
;; or 90% of the steps, whichever is smaller
medchi2 = min(median(chi2[burnndx:nsteps-1,*],dimension=1))

;; the burnndx is the first link where all good chains have been below
;; the median chi2 of the best chain at least once 
burnndxs = lonarr(nchains)
for j=0L, nchains-1 do begin
   tmpndx = (where(chi2[*,j] lt medchi2))(0)
   if tmpndx eq -1 then burnndxs[j] = nsteps-1 $
   else burnndxs[j] = (round(0.1*nsteps) > tmpndx) < round(nsteps*0.9)
endfor

;; only keep chains that increase the total number of links kept
;; after discarding the (increased) burn-in
sorted = sort(burnndxs)
nlinks = lindgen(nchains)
for j=0L, nchains-1 do nlinks[j] = (nsteps-burnndxs[sorted[j]]+1L)*(j+1L)
maxlinks = max(nlinks,ndx)
burnndx = burnndxs[sorted[ndx]]
goodchains = sorted[0:ndx]
goodchains = goodchains[sort(goodchains)]
ngood = ndx+1

;; if less than 3 chains are good, use all chains
;; this isn't going to go well...
if ngood lt 3 then begin
   printandlog, "WARNING: there are insufficient 'good' chains; using them all", logname
   printandlog, "***DO NOT TRUST THESE RESULTS***", logname
   printandlog, "Restart fit at the best-fit found here, refine your starting values, and/or run it much longer", logname
   goodchains = lindgen(nchains)
   ngood = nchains
endif 

;; define the bad chains
if ngood ne nchains then begin
   badchains = sorted[ndx+1:nchains-1] 
   badchains = badchains[sort(badchains)]
endif else badchains = [-1L]

;; never use less than the last 10% of steps
return, burnndx < round(nsteps*0.9)

end
