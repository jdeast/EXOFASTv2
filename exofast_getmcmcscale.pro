;+
; NAME:
;   exofast_getmcmcscale
; PURPOSE:
;   Returns the optimal scale for MCMC chains by varying it until the
;   delta chi^2 = 1.
;
; PROCEDURE:
;   Calculates the chi^2 of the best fit model. Then, for each
;   parameter, takes a small positive step (defined by seedscale) and
;   re-calculates the chi^2. It brackets the value that yields delta
;   chi^2 = 1 by doubling the seedscale until delta chi^2 > 1, then
;   does a binary search to find the value that is as close as
;   possible to delta chi^2 = 1. It repeats with a small negative
;   step. The average of the positive and negative excursion is
;   returned as the optimal scale.
;
;   If the best fit is accurate and the errors are
;   gaussian, uncorrelated, and accurately describe the data
;   (chi^2/dof = 1), this will result in the optimal acceptance rate
;   for MCMC fitting, 44% for 1 parameter, approaching 23% as the
;   number of parameters is large.
;
; INPUTS:
;   BESTPARS   - An NPARS element array of the best fit parameters
;   CHI2FUNC   - A string of the named function that calculates the
;                chi^2 of your model given the parameters specified by
;                BESTPARS.
;
; OPTIONAL INPUTS:
;   SEEDSCALE  - An NPARS arrray that contains the small step that is
;                initially taken. If the chi^2 surface is a smooth,
;                monotonically increasing surface, then it only
;                affects execution time of the program. Since the
;                chi^2 surface is rarely so forgiving, it's better to
;                err on the side of too small. The default is 1d-3 for
;                each parameter.
;   BESTCHI2   - The chi^2 of the best fit parameters. If not given, it
;                is calculated from BESTPARS using CHI2FUNC.
;   ANGULAR    - An array of indices (with up to NPARS elements)
;                indicating which elements are angular values. If
;                specified, and the step sizes determined are larger
;                than 2*pi, the step size will be set to pi.
;                NOTE: angular values must be in radians.
;
; OPTIONAL KEYWORDS:
;  DEBUG       - If set, it will print progress to the screen that
;                will help identify problems.
;
; REVISION HISTORY:
;   2009/11/23 - Written: Jason Eastman - The Ohio State University
;   2013/02/26 - Change print statement to message statement
;   2013/03/08 - If it hits a boundary, take the minimum step as the
;                stepping scale. We just need a good ballpark, DEMC
;                will do the rest.
;   2017/07/19 - Use printandlog for print statements
;-

function exofast_getmcmcscale, bestpars, chi2func, tofit=tofit, $
                               seedscale=seedscale, bestchi2=bestchi2,$
                               angular=angular0, debug=debug, skipiter=skipiter,logname=logname
npars = n_elements(bestpars)
if n_elements(tofit) eq 0 then tofit = indgen(npars)
nfit = n_elements(tofit)
if n_elements(seedscale) ne nfit then seedscale = dblarr(nfit) + 1d-3
maxiter = 1d4

if n_elements(bestchi2) eq 0 then $
  bestchi2 = call_function(chi2func, bestpars)

if ~finite(bestchi2) then begin
   printandlog, 'Best chi^2 is out of bounds; refine bestpars'
   stop
endif

origbestchi2 = bestchi2

;; which values are angular values (error if step > 2*pi)
if n_elements(angular0) eq 0 then angular = [-1] $
else angular = angular0

if keyword_set(debug) then begin
    printandlog, "   Par           Minimum Step                " + $
      "Maximum Step                Value                       " + $
      "Chi^2                    Best Chi^2",logname
endif

mcmcscale = [[seedscale],[seedscale]]

betterfound = 0B

for i=0, nfit-1 do begin

    for j=0,1 do begin

        testpars = bestpars
        
        minstep = 0
        maxstep = 0
        niter = 0
        bestdeltachi2 = !values.d_infinity
        bestscale = 0d0

        repeat begin

            chi2changed = 0

            ;; an infinite step size means it's not constrained
            if ~finite(bestpars[tofit[i]] + mcmcscale[i,j]) or $
               ~finite(bestpars[tofit[i]] - mcmcscale[i,j]) then begin
               printandlog, "EXOFAST_GETMCMCSCALE: Parameter "+strtrim(tofit[i],2) + $
                      " is unconstrained. Check your starting conditions", logname
               chi2 = bestchi2 + 1 
               mcmcscale[i,j] = !values.d_nan
               stop
               goto, next
            endif            
;            if ~finite(bestpars[tofit[i]] + mcmcscale[i,j]) then begin
;            endif

            ;; add the offset to a parameter
            if j eq 0 then begin
                testpars[tofit[i]] = bestpars[tofit[i]] + mcmcscale[i,j]
            endif else begin
                testpars[tofit[i]] = bestpars[tofit[i]] - mcmcscale[i,j]
            endelse

            ;; determine the new chi^2
            chi2 = call_function(chi2func, testpars)
            
            ;; determine the next step size based on deltachi2
            if (chi2 - bestchi2) ge 1.d0 then begin
                ;; if deltachi2 is too large, set max to current value
                maxstep = mcmcscale[i,j]
                ;; the next step will be the midpoint of min/max
                mcmcscale[i,j] = (maxstep + minstep)/2.d0
            endif else if (chi2-bestchi2 ge 0) then begin
                ;; if deltachi2 is too small, set min to current value
                minstep = mcmcscale[i,j]
                ;; if a bound on the max hasn't been determined double the step
                ;; otherwise, take the midpoint of min/max
                if maxstep eq 0 then begin
                    mcmcscale[i,j] *= 2
                endif else mcmcscale[i,j] = (maxstep + minstep)/2.d0
            endif else begin
                if keyword_set(debug) then begin
                    printandlog,'WARNING: better chi2 found by varying parameter '+$
                      strtrim(tofit[i],2) + ' from ' + $ 
                      strtrim(string(bestpars[tofit[i]],format='(f40.10)'),2)+$
                      ' to ' + $
                      strtrim(string(testpars[tofit[i]],format='(f40.10)'),2)+$
                      ' (' + strtrim(chi2,2) + ')',logname
                endif

                if (origbestchi2 - bestchi2) gt 1d0 then betterfound = 1B
                
                ;; chi2 is actually lower! (Didn't find the best fit)
                ;; attempt to fix
                bestpars = testpars

                ;; could be way off, double the step for faster convergence
                mcmcscale[i,j] *= 2d0 
                bestchi2 = chi2
                niter = 0
                chi2changed = 1

            endelse

            deltachi2 = chi2-bestchi2          
            ;; in case we chance upon a better match than we bracket
            ;; (implies chi^2 surface is rough)
            if abs(deltachi2 - 1d0) lt abs(bestdeltachi2 - 1d0) then begin
                bestdeltachi2 = deltachi2
                bestscale = mcmcscale[i,j]
            endif

            ;; can't always sample fine enough to get exactly
            ;; deltachi2 = 1 because chi^2 surface not perfectly smooth
            if abs(minstep - maxstep) lt 1d-12 or niter gt maxiter then begin
               testpars[tofit[i]] = bestpars[tofit[i]] - 2d0*mcmcscale[i,j]
               lowboundchi2 = call_function(chi2func, testpars)
               testpars[tofit[i]] = bestpars[tofit[i]] + 2d0*mcmcscale[i,j]
               hiboundchi2 = call_function(chi2func, testpars)

               if ~finite(chi2) or ~finite(lowboundchi2) or ~finite(hiboundchi2) then begin
                  if j eq 0 then bound = 'upper' $
                  else bound = 'lower'
;                  printandlog, 'Reached the ' + bound + ' bound before hitting deltachi^2 = 1 for  parameter ' + strtrim(tofit[i],2) + '.'
;                  printanglog, 'Using a scale of 1/10 the distance to the boundary (' + strtrim(bestscale/10d0,2) + ').'
                  printandlog, 'The ' + bound + ' bound for parameter ' + strtrim(tofit[i],2) + ' is critical; it must be physically and independently motivated.'
                  if bestscale ne 0d0 then $
                     mcmcscale[i,j] = bestscale/100d0 $
                  else mcmcscale[i,j] = maxstep/100d0;!values.d_nan
                  chi2 = bestchi2 + 1 
                  goto, next
               endif else if not chi2changed then begin
                  if abs(bestdeltachi2 - 1.d0) lt 0.75 then begin
                     mcmcscale[i,j] = bestscale
                  endif else begin
                     if bestdeltachi2 eq 0d0 then newscale = bestscale/100d0 $
                     else newscale = bestscale/bestdeltachi2/10d0 ;; better to err on the side of too small
                     printandlog, 'Cannot find the value for which deltachi^2 = 1 for parameter ' +$
                                  strtrim(tofit[i],2) + '; assuming rough chi^2 surface. Using delta chi^2 = ' +$
                                  strtrim(bestdeltachi2,2) + ' and scaling step (' + strtrim(bestscale,2) +$
                                  ') to ' + strtrim(newscale,2),logname
                     
                     ;; extrapolate scale to delta chi^2 = 1 
                     ;; (may cause issues near boundaries)
                     mcmcscale[i,j] = newscale
                  endelse
               endif
               chi2 = bestchi2 + 1d0
            endif

            ;; if the parameter has no influence on the chi2
            if abs(chi2 - bestchi2 - 1.d0) eq 1.d0 and niter gt maxiter then $
               printandlog, 'ERROR: changing parameter ' + strtrim(tofit[i],2) + $
                            ' does not change the chi^2. Exclude from fit.', logname

            ;; if angle is so poorly constrained 
            ;; no value has delta chi^2 = 1
            if (where(angular eq tofit[i]))(0) ne -1 and $
              minstep gt (2.d0*!dpi) then begin 
               printandlog, 'WARNING: no constraint on angle',logname
                mcmcscale[i,j] =  !dpi
                chi2 = bestchi2 + 1.d0
            endif

            ;; near a boundary
            ;; exclude this direction from the scale calculation
;            if ~finite(chi2) then begin
;                chi2 = bestchi2 + 1.d0
;                mcmcscale[i,j] = minstep
;            endif
            niter++

            ;; print the steps and chi^2 it's getting
            if keyword_set(debug) then begin
                if j eq 0 then str = '(hi)' $
                else str = '(lo)'                
                printandlog, string(tofit[i], str, minstep, maxstep, $
                  testpars[tofit[i]],chi2,bestchi2, $
                  format='(i3,x,a4,x,f0.18,x,f0.18,x,f0.18,x,f0.18,x,f0.18)'),logname
             endif
            next:

         endrep until abs(chi2 - bestchi2 - 1.d0) lt 1d-8
    endfor    
endfor

;; let's iterate to see if we can do better
if betterfound and ~keyword_set(skipiter) then begin
   printandlog, 'Better Chi^2 found (' + strtrim(bestchi2,2) + ' vs ' + strtrim(origbestchi2,2) + '. Restarting EXOFAST_GETMCMCSCALE at new best fit',logname
   printandlog, 'This means the best fit was not properly identified and MCMC results may be slow or suspect!',logname
   printandlog, 'If the chi^2 difference is large, you may want to revise your starting parameters and/or check your priors, then restart',logname
   printandlog, 'Proceed with caution and skepticism',logname
   return, exofast_getmcmcscale(bestpars, chi2func, tofit=tofit, $
                                seedscale=seedscale, bestchi2=bestchi2,$
                                angular=angular0, debug=debug, skipiter=skipiter)
endif

;; replace undefined errors with the other
bad = where(~finite(mcmcscale[*,0]) or mcmcscale[*,0] eq 0d0)
if bad[0] ne -1 then mcmcscale[bad,0] = mcmcscale[bad,1]
bad = where(~finite(mcmcscale[*,1]) or mcmcscale[*,1] eq 0d0)
if bad[0] ne -1 then mcmcscale[bad,1] = mcmcscale[bad,0]
bad = where(~finite(mcmcscale) or mcmcscale eq 0d0,nbad)

if bad[0] ne -1 then stop;return, -1

;; return the average of high and low
return, total(mcmcscale,2)/2.d0

end
