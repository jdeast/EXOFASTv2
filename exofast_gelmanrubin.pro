;+
; NAME:
;   exofast_gelmanrubin
;
; PURPOSE: 
;   Calculates the Gelman-Rubin Statistic for several parallel Markov
;   chains, and determines whether or not the chains are well mixed
;   according to the Ford 2006 prescription: 
;   http://adsabs.harvard.edu/abs/2006ApJ...642..505F 
;
;   This is meant to be a dependency of EXOFAST_DEMC. An example of
;   its use can be seen there.
;
; PROCEDURE
;   Calculates the Gelman-Rubin statistic and the number of
;   independent draws for each parameter, as defined by Ford, 2006. The
;   chain is considered well-mixed if all parameters have a
;   Gelman-Rubin statistic of <= MAXGR (default = 1.01) and >= MINTZ
;   (default = 1000) independent draws.
;
; CALLING SEQUENCE:
;   ismixed = exofast_gelmanrubin(pars)
;
; INPUTS:
;   PARS0       - A 3 dimensional array (NPARS,NSTEPS,NCHAINS) of
;                 parameter values
; OPTIONAL INPUTS:
;   MAXGR       - The maximum Gelman Rubin statistic that is
;                 considered well-mixed (default=1.01)
;   MINTZ       - The minimum number of independent draws that is
;                 considered well-mixed (default=1000)
;
; OPTIONAL OUTPUTS:
;   GELMANRUBIN - An NPARS element array containing the Gelman-Rubin
;                 statistic for each parameter (equation 25).
;   TZ          - An NPARS element array containing the number of
;                 independent draws for each parameter (equation 26).
;
; OPTIONAL KEYWORDS:
;   ANGULAR     - An array of indices (with up to NPARS elements)
;                 indicating which elements are angular values. If
;                 specified, these parameters will be re-centered
;                 about the mode value before calculating the
;                 statisitics (the original array will not be
;                 modified). If left blank, no scaling will be performed. 
;                 NOTE: angular values must be in radians.
;  
; REVISION HISTORY:
;   2010/03/01 - Written: Jason Eastman - The Ohio State University
;   2013/03/18 - Added MAXGR, MINTZ optional inputs.
;-

function exofast_gelmanrubin, pars0, gelmanrubin, tz, angular=angular, maxgr=maxgr, mintz=mintz

if n_elements(mintz) eq 0 then mintz = 1000
if n_elements(maxgr) eq 0 then maxgr = 1.01

pars = pars0 ;; don't modify input parameters

sz = size(pars)
if sz[0] ne 3 then message, 'ERROR: pars must have 3 dimensions'

npars = double(sz[1])
nsteps = double(sz[2])
nchains = double(sz[3])

if nsteps eq 1 then message, 'ERROR: NSTEPS must be greater than 1'
if nchains lt 3 then message, 'ERROR: NCHAINS must be greater than 2'

;; If an angular parameter, must be careful calculating the median.
;; Center about the mode first
nangular = n_elements(angular)
if nangular ne 0 then begin
    if angular[0] ne -1 then begin
        for i=0, nangular-1 do begin
            hist = histogram(pars[angular[i],*,*],nbins=100,locations=x)
            max = max(hist,modendx)
            mode = x[modendx]
            
            for j=0, nchains-1 do begin
                toohigh = where(pars[angular[i],*,j] gt (mode + !pi))
                if toohigh[0] ne -1 then pars[angular[i],toohigh,j] -= 2.d0*!dpi
                
                toolow = where(pars[angular[i],*,j] lt (mode - !pi))
                if toolow[0] ne -1 then pars[angular[i],toolow,j] += 2.d0*!dpi
            endfor
        endfor
    endif
endif

;; Equation 21: W(z) in Ford 2006
variances = dblarr(npars,nchains)
for i=0,npars-1 do for k=0, nchains-1 do variances[i,k]=stddev(pars[i,*,k])^2
meanofvariances = total(variances,2)/nchains

;; Equation 23: B(z) in Ford 2006
means = total(pars,2)/nsteps
varianceofmeans = total((means-replicate(1,nchains)##$
                         total(means,2)/nchains)^2,2)/(nchains-1)

bz = varianceofmeans*nsteps

;; Equation 24: varhat+(z) in Ford 2006
varz = (nsteps-1.d0)/nsteps*meanofvariances + varianceofmeans

;; Equation 25: Rhat(z) in Ford 2006
gelmanrubin = sqrt(varz/meanofvariances)

;; Equation 26: T(z) in Ford 2006
tz = nchains*nsteps*(varz/bz < 1)

;; well-mixed criteria
if min(tz,/nan) gt mintz and max(gelmanrubin,/nan) lt maxgr then return, 1

;; not well mixed
return, 0

;; the slow (but more comprehensible) way...
zstarc = dblarr(npars,nchains)
wz = dblarr(npars)
zstarstar = dblarr(npars)
bz = dblarr(npars)
varz2 = dblarr(npars)
rz = dblarr(npars)

for i=0L, npars-1 do begin

    for j=0L, nsteps-1 do begin
        for k=0L, nchains-1 do begin
            zstarc[i,k] += pars[i,j,k]                   ;; equation 20
            zstarstar[i] += pars[i,j,k]                  ;; equation 22
        endfor
    endfor
    zstarc[i,*] /= nsteps
    zstarstar[i] /= (nsteps*nchains)

    for j=0L, nsteps-1 do begin
        for k=0L, nchains-1 do begin
            wz[i] += (pars[i,j,k] - zstarc[i,k])^2       ;; equation 21
        endfor
    endfor
    wz[i] /= (nchains*(nsteps-1.d0))

    for k=0L, nchains-1 do begin
        bz[i] += (zstarc[i,k] - zstarstar[i])^2          ;; equation 23        
    endfor
    bz[i] = bz[i]*nsteps/(nchains-1.d0)

    varz2[i] = wz[i]/nsteps*(nsteps-1.d0) + bz[i]/nsteps ;; equation 24
    rz[i] = sqrt(varz2[i]/wz[i])                         ;; equation 25
endfor


end
