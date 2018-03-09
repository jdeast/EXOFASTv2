;+
; NAME:
;   exofast_ptmc
; PURPOSE:
;   A Parallel Tempering, Differential Evolution, Markov Chain Monte
;   Carlo implementation to fit data and return the parameter
;   distributions, robust to widely-separated local minima.
;
; DESCRIPTION:
;   Begins by determining the correct stepping scale by finding
;   differences in each parameters that yield a delta chi^2 = 1, then
;   uses 5*scale*randomn offsets in each parameter for the starting
;   points of each chain.
;
;   It then begins a Differential Evolution Markov Chain Monte Carlo
;   fit (ter Braak, 2006)
;   http://www.stat.columbia.edu/~gelman/stuff_for_blog/cajo.pdf. The
;   only slight modification to this basic alogorithm is that we
;   determine the magnitude of the uniform random deviate for each
;   parameter. We have found the dynamic range of our parameters can
;   be large, and thus not well suited to the one-size fits all
;   approach recommended there.
;
;   It does this for each of NTEMPS "temperatures", ranging from 1 to
;   TF, each of which is successively more likely to accept a worse
;   step. Swaps between simulations at adjacent temperatures are
;   accepted with the Metropolis-Hastings algorithm. This allows the
;   algorithm to sample widely-separated local minima
;   simultaneously. The parameters returned only come from the T=1
;   chain. By default, NTEMPS=1, and no parallel tempering is done
;   (for well-behaved parameter distributions, this is a waste of
;   time). The hot parameters can optionally be returned in the
;   HOTPARS array, but they should not be used to infer the posterior
;   and it substantially increases the memory footprint of the
;   program.
;
;   After taking 5% of MAXSTEPS, the program will estimate how many steps
;   will be necessary to be well-mixed. If it is not expected to be
;   well-mixed before taking the maximum number of steps, it will
;   output a warning and a recommended setting for NTHIN. This should
;   be accurate at the factor of 2-3 level.
;
;   The program stops when the chains are well-mixed, as defined by
;   Ford 2006 (http://adsabs.harvard.edu/abs/2006ApJ...642..505F)
;   using the Gelman-Rubin Statistic and number of independent draws,
;   or when each chain has taken MAXSTEPS, whichever is first.
;
;   Every step in the chain before all chains have crossed below the
;   median chi^2 will be considered the (the "burn-in"), and not used
;   for convergence tests. BURNNDX will specify the index of the first
;   useable chain so as not to be biased by the starting criteria.
;
; CALLING SEQUENCE:
;   exofast_demc, bestpars, 'mychi2', pars [,CHI2=, TOFIT=,$
;                 SCALE=, SEED=,RANDOMFUNC=,NTHIN=, MAXSTEPS=,/DONTSTOP,$
;                 NCHAINS=, ANGULAR=,BURNNDX=,/REMOVEBURN]
;
; INPUTS:
;   BESTPARS   - Array of the best-fit parameters. An accurate initial
;                fit is required to find the correct step size. MCMC
;                is not ideal for fitting data, just characterizing
;                the errors.
;   CHI2FUNC   - A string that specifies the name of the user-defined
;                function that calculates the chi^2 of a given
;                parameter set. Required data should be put in a
;                COMMON block. Its calling sequence must be:
;                chi2 = chi2func(pars, determinant=determinant)
;
; OPTIONAL INPUTS:
;   TOFIT      - An array that indexes the parameters to fit. If not
;                specified, all parameters are fit.
;   ANGULAR    - An array that indexes the parameters that are
;                angular (must be in radians). This will enable
;                special ways to calculate the median and convergence
;                statistics, which may otherwise fail. Default is none.
;   SCALE      - An NPARS array containing the stepping scales for
;                each parameter. This is not recommended for normal
;                use. If not specified, it will be automatically
;                determined using EXOFAST_GETMCMCSCALE.
;   NTHIN      - Saves only every NTHINth link in each chain. This
;                is only recommended when the autocorrelation between
;                steps is large and memory management issues may
;                arise. Results will be more accurate the smaller this
;                is -- as long as the chain is well-mixed (pay
;                attention to the warnings).
;   MAXSTEPS   - The maximum number of steps (after thinning) for each
;                chain to take.  The default is 100,000. Take care when
;                increasing this number, as memory allocation issues
;                may arise (especially with 32 bit machines). If the
;                chain is not well-mixed, it is usually better to
;                increase NTHIN instead.
;   NCHAINS    - The number of independent chains to run. The
;                execution time scales ~linearly with this number, and
;                must be at least 3 to calculated convergence
;                statistics. The default is 10. 
;   SEED       - A seed for the random number generator. Do not mix
;                seeds from different random number
;                generators, and do not use several different seeds
;                for the same generator. If you use a random number
;                generator elsewhere, you should call EXOFAST_DEMC with
;                SEED. Default is -systime(/seconds). 
;   RANDOMFUNC - A string for the name of the random number generator
;                function. You can use your own as long as it returns
;                a uniform double precision deviate between 0 and 1,
;                accepts the /NORMAL keyword to generate random
;                gaussian deviates, and can return at least a two
;                dimensional array of deviates. You may also choose
;                between "EXOFAST_RANDOM" (default) or "RANDOMU".
;              - EXOFAST_RANDOM is an IDL wrapper for pg_ran. This is
;                a state of the art, bullet-proof generator with no
;                known statistical flaws, based on the third edition
;                of numerical recipes, pg 342. It is ~120x slower than
;                RANDOMU, though for typical DEMC fits, this accounts
;                for ~1% of the total runtime.
;              - RANDOMU is IDL's built-in, ran1-based generator with a
;                periodicity of ~10^9. ran1 is formally admonished by
;                its author, Numerical Recipes, but shown to give
;                identical results as the state of the art generator
;                See Eastman et. al., 2013 for more details.
;                http://adsabs.harvard.edu/abs/2013PASP..125...83E
;   MAXGR      - The maximum Gelman Rubin statistic that is
;                considered well-mixed (default=1.01)
;   MINTZ      - The minimum number of independent draws that is
;                considered well-mixed (default=1000)
;   NTEMPS     - The number of chains to run at different
;                temperatures. The more temperatures, the
;                easier it is to make swaps between temperatures, but
;                the longer it takes. The default is 1 (no parallel
;                tempering). 8 is a good value. Note that only the T=1
;                chain is kept.
;   TF         - The temperature of the hottest chain. The higher the
;                temperature, the more parameter space it will probe,
;                but the chains will be less likely to swap. Default is
;                200, which can robustly find minima separated by 50
;                sigma. A value of 1 will result in no temperature
;                difference between chains.
;   HOTPARS    - By default, the hot chains from the parallel
;                tempering are dynamically discarded
;                for efficiency. If this parameter is set, the
;                parameters from the hot chains will be stored
;                here. This is not useful to determine the posterior,
;                and should only be used for debugging. Note this
;                increases the (already-substantial) memory footprint
;                of the program by roughly a factor of NTEMPS.
;
; OPTIONAL KEYWORDS:
;   DONTSTOP   - If set, the program will not stop when it is
;                converged, but will instead take MAXSTEPS steps in
;                each chain.
;   REMOVEBURN - If set, the burn-in will be removed in the returned
;                parameter array.
; OUTPUTS:
;   PARS       - NFIT x NSTEPS x NCHAINS array of parameters.
;                NOTE: NSTEPS will be the lesser of when it is
;                converged or MAXSTEPS. 
;   CHI2       - NSTEPS x NCHAINS array with the chi^2 for each set of
;                parameters
;   DERIVED    - An NDERIVED x NSTEPS x NCHAINS array of parameters
;                that are derived during the chi^2 function. To
;                conserve memory, these should be reserved for
;                parameters that cannot be trivially rederived.
;
; REVISION HISTORY:
;   2012/06 - Public Release - Jason Eastman (LCOGT)
;   2012/12 - When parameters aren't mixed, display the index
;             to the (more) static PARS array, not the index to the
;             dynamic TOFIT array.
;           - Removed memory-hungry vestigial code
;   2013/03 - Add GELMANRUBIN and TZ outputs, MAXGR, MINTZ inputs.
;             Added DERIVED output, to return non-trivially derived
;             parameters.
;             Added check for errant chains (stuck in local minima)
;             and discard them if the chains are not well-mixed at the
;             end. This is sketchy, but better than nothing.
;             This was a relatively common failure mode.
;   2013/11 - MAJOR CHANGE -- added parallel tempering option
;-

pro exofast_ptmc, bestpars,chi2func,pars,chi2=chi2, tofit=tofit,$
                  scale=scale,seed=seed,randomfunc=randomfunc,$
                  nthin=nthin, maxsteps=maxsteps, dontstop=dontstop,$
                  nchains=nchains, angular=angular,$
                  burnndx=burnndx, removeburn=removeburn, $
                  gelmanrubin=gelmanrubin, tz=tz, maxgr=maxgr, mintz=mintz, derived=derived, $
                  ntemps=ntemps, tf=tf, hotpars=hotpars

;; default values
if n_elements(maxsteps) eq 0 then maxsteps = 1d5
if n_elements(nthin) eq 0 then nthin = 1L
if n_elements(tofit) eq 0 then tofit = indgen(n_elements(bestpars))
if n_elements(mintz) eq 0 then mintz = 1000d0
if n_elements(maxgr) eq 0 then maxgr = 1.01d0
nfit = n_elements(tofit)
!except = 0 ;; don't display errors for NaN and infinity

;; the number of parallel temperatures to run
if n_elements(ntemps) eq 0 then begin
   temps = 1d0
   ntemps = 1d0
endif else begin
   if n_elements(tf) eq 0 then tf = 200d0
   temps = (1d0/tf)^(dindgen(ntemps)/(ntemps-1d0))
endelse

;; defaults for Differential Evolution MCMC
gamma = 2.38d0/sqrt(2d0*nfit)
if n_elements(nchains) eq 0 then nchains=2d0*nfit

;; use robust generator if estimated steps is large 
;; NOTE: necessarily ignores any random numbers used by chi^2 routine
if n_elements(randomfunc) eq 0 then begin
   nfit = n_elements(tofit)
   ;; number of random deviates required
   ;; Gaussian deviates are made in pairs, and 4/pi must be made for each
   nrand = ((ceil(nfit/2d0)*2d0)*(4d0/!dpi)+1d0)*maxsteps*nchains*nthin
   if nrand le 5d7 then randomfunc = 'randomu' $
   else randomfunc = 'exofast_random'
endif

;; angular must be indexed within tofit to use with gelman-rubin
nang = n_elements(angular)
if nang gt 0 then gelmanangular=intarr(nang)
for i=0, nang-1 do gelmanangular[i] = (where(tofit eq angular[i]))(0)

;; initialize arrays
pars = dblarr(nfit,maxsteps,nchains,ntemps)
chi2 = dblarr(maxsteps,nchains,ntemps)
niter = 1L
newpars = bestpars
oldpars = replicate(1d0,ntemps)##bestpars
olddet = dblarr(nchains,ntemps)

;; get the MCMC step scale for each parameter
if n_elements(scale) eq 0 then $
   scale = exofast_getmcmcscale(bestpars,chi2func,tofit=tofit,angular=angular)
if scale[0] eq -1 then begin
   message, 'No scale found',/continue
   pars = -1
   return
endif

;; initialize each chain
for j=0, nchains-1 do begin
   for k=0, ntemps-1 do begin
      ;; repeat until a finite chi^2 
      ;; i.e., don't start outside a boundary; you'll never get back
      niter = 0d0
      repeat begin
         if j eq 0 then factor = 0d0 else factor = (sqrt(500d0/nfit) > 3d0) ;; first chain starts at best fit
         ;; start 5 steps from best value
         pars[0:nfit-1,0,j,k] = bestpars[tofit] + $
                                factor/exp(niter/1000d)*scale*call_function(randomfunc,seed,nfit,/normal)
         newpars[tofit] = pars[0:nfit-1,0,j,k]

         ;; find the chi^2
         chi2[0,j,k] = call_function(chi2func, newpars, determinant=det, derived=dpar)
         pars[*,0,j,k]=newpars[tofit] ;; in case chi2 function changes pars
         olddet[j,k] = det
         niter += 1d0
      endrep until finite(chi2[0,j,k])
      if n_elements(dpar) ne 0 then begin
         if j eq 0 then begin
            nderived = n_elements(dpar)
            derived = dblarr(nderived,maxsteps,nchains,ntemps)
         endif
         derived[*,0,j,k] = dpar
      endif
   endfor
endfor
if n_elements(dpar) ne 0 then oldderived = reform(derived[*,0,*,*],nderived,nchains,ntemps)

nextrecalc = 100L
npass = 1L 
nstop = 0L
naccept = 1d0
tz0 = 0d0
tzsteps = 0L
alreadywarned = 0L
t0 = systime(/seconds)
nswap = 0L
nswapattempt = 0L

;; start MCMC chain
for i=1L,maxsteps-1L do begin
   for j=0L, nchains-1L do begin

      oldpars[tofit,*] = pars[*,i-1,j,*]
      oldchi2 = chi2[i-1,j,*]
      
      ;; automatically thin the chain (saves memory)
      for k=0L, nthin-1L do begin

         for m=0L, ntemps-1L do begin
            
            ;; differential evolution mcmc step -- see Ter Braak, 2006
            repeat r1=floor(call_function(randomfunc,seed)*(nchains)) $
            until r1 ne j ;; a random chain
            repeat r2=floor(call_function(randomfunc,seed)*(nchains)) $
            until r2 ne j and r2 ne r1 ;; another random chain
            newpars[tofit] = oldpars[tofit,m] + $ 
                             gamma*(pars[*,i-1,r1,m]-pars[*,i-1,r2,m]) + $
                             (call_function(randomfunc,seed,nfit)-0.5d0)*scale/10d0
         
            ;; calculate the chi^2 of the new step
            newchi2 = call_function(chi2func,newpars,determinant=det,derived=dpar)
            C = olddet[j,m]/det*exp(temps[m]*(oldchi2[m] - newchi2)/2d0)
            
            ;; accept the step; update values
            if call_function(randomfunc,seed) lt C then begin
               if m eq 0 then naccept++
               olddet[j,m] = det
               oldpars[*,m] = newpars
               oldchi2[m] = newchi2
               if n_elements(dpar) ne 0 then oldderived[*,j,m] = dpar
            endif ;; else keep previous values

            ;; Parallel Tempering
            ;; http://www.lindonslog.com/programming/openmp/parallel-tempering-algorithm-c/
            ;; attempt swap between temperatures with probability 1/2
            if call_function(randomfunc,seed) lt 0.5d0 and m ne 0 then begin

               nswapattempt++
               ;; accept the swap with Metropolis-Hastings criteria
               C = exp((oldchi2[m-1]-newchi2)/2d0) ;; do I need a determinant??
               if call_function(randomfunc,seed) lt C then begin
                  tmp = oldpars[*,m-1]
                  oldpars[*,m-1] = newpars
                  oldpars[*,m] = tmp

                  ;; don't forget to swap auxiliary arrays, too
                  tmp = oldchi2[m-1]
                  oldchi2[m-1] = oldchi2[m]
                  oldchi2[m] = tmp
                 
                  tmp = olddet[j,m-1]
                  olddet[j,m-1] = olddet[j,m]
                  olddet[j,m] = tmp

                  if n_elements(dpar) ne 0 then begin
                     tmp = oldderived[*,j,m-1]
                     oldderived[*,j,m-1] = oldderived[*,j,m]
                     oldderived[*,j,m] = tmp
                  endif
                  nswap++
               endif
            endif

         endfor
      
         ;; update the link in the chain
         pars[*,i,j,*] = oldpars[tofit,*]
         chi2[i,j,*] = oldchi2
         if n_elements(dpar) ne 0 then derived[*,i,j,*] = oldderived[*,j,*]

      endfor
   endfor
   
   ;; Test for convergence as outlined in Ford 2006
   ;; must be converged for 6 consecutive passes
   ;; tz > 1000 and Gelman-Rubin < 1.01 => converged 
   if i eq nextrecalc then begin
      ;; discard the burn-in (the first point at which all chains
      ;; have had at least one chi^2 lower than the median chi^2)

      ;; the minimum chi2 of each chain
      minchi2 = min(chi2[0:i,*,0],dimension=1)

      ;; the median chi2 of the best chain
      medchi2 = min(median(chi2[0.1*i:i,*,0],dimension=1))

      ;; good chains must have some values below this median
      goodchains = where(minchi2 lt medchi2,ngood)
      if ngood lt 3 then begin
         goodchains = lindgen(nchains)
         ngood = nchains
      endif

      burnndx = 0L
      for j=0L, ngood-1 do begin
         tmpndx = (where(chi2[0:i,goodchains[j],0] lt medchi2))(0)
         if tmpndx gt burnndx then burnndx = tmpndx
      endfor
      ;; allows Gelman-Rubin calculation if one chain is being problematic
      burnndx = burnndx < (i-3) 
      
      ;; calculate the Gelman-Rubin statistic (remove burn-in)
      converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:i,goodchains,0],$
                                      gelmanrubin,tz,angular=gelmanangular,$
                                      mintz=mintz,maxgr=maxgr)
      
      ;; estimate the number of steps it will take until it is converged
      tz0 = [tz0,min(tz)]
      tzsteps = [tzsteps,i]
      ntz = n_elements(tz0)
      if not alreadywarned then begin
         if i gt maxsteps/20d0 and ntz gt 3 then begin
            ;; fit a line to the number of independent draws
            coeffs = poly_fit(tzsteps[1:ntz-1],tz0[1:ntz-1],1,yfit=fit)
            
            ;; extrapolate to 2000, when it'll be considered converged 
            ;; MINTZ is considered converged, but this is far from exact...
            stepstoconvergence = (2d0*mintz - coeffs[0])/coeffs[1]
            
            ;; if it won't be converged, warn the user
            if stepstoconvergence gt maxsteps then begin
               bestnthin = round(stepstoconvergence*nthin/maxsteps)
               print,'WARNING: The chain is not expected to be well-mixed; '+$
                     'set NTHIN to >~ ' + strtrim(bestnthin,2) + ' to fix'
            endif else if stepstoconvergence le 0 then begin
               print,'WARNING: The chains are not mixing. This may be a while.'
            endif else print, string(stepstoconvergence*50d0/maxsteps,$
                     format='("EXOFAST_DEMC: The chain is expected to be ' + $
                              'well-mixed after ",i2,"% completed     ")')
            alreadywarned = 1 ;; only display message once
         endif
      endif
      
      ;; must pass 6 consecutive times
      if converged then begin
         if nstop eq 0 then nstop = i
         nextrecalc = long(nstop/(1.d0-npass/100.d0))
         npass++
         if npass eq 6 then begin
            if not keyword_set(dontstop) then break
            nextrecalc = maxsteps
         endif
      endif else begin
         nextrecalc = long(i/0.9d0)
         nstop = 0L
         npass = 1L
      endelse
   endif
   
   ;; print the progress, cumulative acceptance rate, and time remaining
   acceptancerate = strtrim(string(naccept/double(i*nchains*nthin)*100,$
                                   format='(f6.2)'),2)
   swaprate = strtrim(string(nswap/double(nswapattempt)*100,$
                             format='(f6.2)'),2)
   timeleft = (systime(/seconds)-t0)*(maxsteps/(i+1d)-1d)
   units = ' seconds '
   if timeleft gt 60 then begin
      timeleft /= 60
      units = ' minutes '
      if timeleft gt 60 then begin
         timeleft /= 60
         units = ' hours   '
         if timeleft gt 24 then begin
            timeleft /= 24
            units = ' days    '
         endif
      endif
   endif

   ;; Windows formatting is messy, only output every 1%
   if !version.os_family ne 'Windows' or $
      i mod round(maxsteps/1000d0) eq 0 then begin
      timeleft = strtrim(string(timeleft,format='(f255.2)'),2)
      if finite(swaprate) then begin
         format='("EXOFAST_DEMC:",f6.2,"% done; acceptance rate=",a,"%; swap rate=",a,"%; ' + $ 
                'time left: ",a,a,$,%"\r")'
         print, format=format,100.d0*(i+1)/maxsteps,acceptancerate,swaprate,timeleft,units
      endif else begin
         ;; no parallel tempering, don't print swap information
         format='("EXOFAST_DEMC:",f6.2,"% done; acceptance rate=",a,"%; ' + $ 
                'time left: ",a,a,$,%"\r")'
         print, format=format,100.d0*(i+1)/maxsteps,acceptancerate,timeleft,units
       endelse
   endif

endfor
print ;; don't overwrite the final line
print ;; now just add a space

;; if didn't converge, or user didn't want to stop use all steps
if npass ne 6 or keyword_set(dontstop) then nstop = maxsteps-1

;; the minimum chi2 of each chain
minchi2 = min(chi2[0:nstop,*,0],dimension=1)

;; the median chi2 of the best chain
medchi2 = min(median(chi2[0.1*nstop:nstop,*,0],dimension=1))

;; good chains must have some values below this median
goodchains = where(minchi2 lt medchi2,ngood)

burnndx = 0L
for j=0L, ngood-1 do begin
   tmpndx = (where(chi2[0:nstop,goodchains[j],0] lt medchi2))(0)
   if tmpndx gt burnndx then burnndx = tmpndx
endfor
;; allows Gelman-Rubin calculation if one chain is being problematic
burnndx = burnndx < (i-3) 

;; if there are bad chains, is it better without them?
;; (it usually is unless it's a very short run)
if ngood ne nchains or npass ne 6 then begin
   ;; are all the chains mixed?
   converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,*,0],$
                                   gelmanrubin,tz,angular=gelmanangular, $
                                   maxgr=maxgr, mintz=mintz)
   
   if ngood gt 2 then begin
      converged2 = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,goodchains,0],$
                                       gelmanrubin2,tz2,angular=gelmanangular, $
                                       maxgr=maxgr, mintz=mintz)
      bad = where(tz lt mintz or gelmanrubin gt maxgr)
      if bad[0] ne -1 then extrabad = where(tz2[bad] lt tz[bad] or $
                            gelmanrubin2[bad] gt gelmanrubin[bad]) $
      else extrabad = -1

   endif else begin
      bad = where(tz lt mintz or gelmanrubin gt maxgr)
      extrabad = -1
   endelse

   if bad[0] eq -1 then begin
      if npass ne 6 then print, 'Chains are marinally well-mixed' $
      else print,'All chains are mixed, but some chains never got below '+$
            'the median chi^2. Is this even possible?'
   endif else if extrabad[0] eq -1 and ngood ne nchains and ngood gt 2 then begin
      ;; without errant chains, it's better
      print, 'Discarding ' + strtrim(round(nchains-ngood),2) + '/'+$
             strtrim(round(nchains),2) + ' chains stuck in local minima.'
      pars = pars[*,*,goodchains,*]
      chi2 = chi2[*,goodchains,*]
      if n_elements(derived) ne 0 then derived = derived[*,*,goodchains]
      nchains = ngood
      tz = tz2
      gelmanrubin = gelmanrubin2
      bad = where(tz lt mintz or gelmanrubin gt maxgr)
      if not converged2 then begin
         print, 'Remaining chains are still not mixed. Must run longer.'
         print, '       Parameter  Gelman-Rubin Independent Draws'
         forprint, tofit[bad], gelmanrubin[bad],tz[bad],textout=1
      endif
   endif else begin
      print, 'No obviously-errant chains. Must run longer.'
      print, '       Parameter  Gelman-Rubin Independent Draws'
      forprint, tofit[bad], gelmanrubin[bad],tz[bad],textout=1
   endelse
endif


if 0 then begin
;; if hit maxsteps and never converged
if npass ne 6 then begin

   converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,goodchains,0],$
                                   gelmanrubin,tz,angular=gelmanangular, $
                                   maxgr=maxgr, mintz=mintz)
   bad = where(tz lt mintz or gelmanrubin gt maxgr)

   if bad[0] ne -1 then begin
      message, 'WARNING: The Gelman-Rubin statistic indicates ' + $
               'the following parameters are not well-mixed', /continue
      print, '       Parameter  Gelman-Rubin Independent Draws'
      forprint, tofit[bad], gelmanrubin[bad],tz[bad],textout=1
   endif else message, 'WARNING: The chain did not pass 6 consecutive ' + $
                       'tests and may be marginally well-mixed.',/continue

   ;; Maybe it's just an errant chain or two. See if discarding them helps. 
   minchi2 = min(chi2[0:nstop,*,0],dimension=1)

   ;; the median chi2 of the best chain
   medchi2 = min(median(chi2[0.1*nstop:nstop,*,0],dimension=1))
   
   ;; good chains must have some values below this median
   goodchains = where(minchi2 lt medchi2,ngood)

   if ngood ge 3 and ngood lt nchains then begin
      ;; recompute the burn index
      medchi2 = median(chi2[0:nstop,goodchains,0])
      burnndx2 = 0L
      for j=0L, ngood-1 do begin
         tmpndx = (where(chi2[0:nstop,goodchains[j],0] lt medchi2))(0)
         if tmpndx gt burnndx2 then burnndx2 = tmpndx
      endfor
      burnndx2 = burnndx2 < (maxsteps-3) 

      ;; now are they converged?
      converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,goodchains,0],$
                                      gelmanrubin2,tz2,angular=gelmanangular, $
                                      maxgr=maxgr, mintz=mintz)
      bad = where(tz2 lt mintz or gelmanrubin2 gt maxgr)
      if bad[0] ne -1 then begin
         extrabad = where(tz2[bad] lt tz[bad] or $
                          gelmanrubin2[bad] gt gelmanrubin[bad]) 
      endif else extrabad = -1

      if extrabad[0] eq -1 and ngood ne nchains then begin
         ;; without errant chains, it's better
         print, 'Discarding ' + strtrim(round(nchains-ngood),2) + '/'+$
                strtrim(round(nchains),2) + ' chains stuck in local minima.'
         pars = pars[*,*,goodchains,*]
         chi2 = chi2[*,goodchains,*]
         if n_elements(derived) ne 0 then derived = derived[*,*,goodchains]
         nchains = ngood
         tz = tz2
         gelmanrubin = gelmanrubin2
         burnndx = burnndx2
         if bad[0] eq -1 then begin
            print,'Remaining chains are mixed, but interpret results with care.'
         endif else begin
            print, 'Remaining chains are still not mixed.'
            print, '       Parameter  Gelman-Rubin Independent Draws'
            forprint, tofit[bad], gelmanrubin[bad],tz[bad],textout=1
         endelse     
      endif else print, 'No obviously-errant chains. Must run longer.'
   endif else print, 'No obviously-errant chains. Must run longer.'

endif
endif


;; remove the burn-in/uncalculated parameters
;; keep them by default so people can do their own analysis
if keyword_set(removeburn) then begin
   pars = pars[*,burnndx:nstop,*,*]
   chi2 = chi2[burnndx:nstop,*,*]
   if n_elements(derived) ne 0 then derived = derived[*,burnndx:nstop,*,*]
endif else begin
   pars = pars[*,0:nstop,*,*]
   chi2 = chi2[0:nstop,*,*]
   if n_elements(derived) ne 0 then derived = derived[*,0:nstop,*,*]
endelse

;; calculate the runtime and units
runtime = systime(/seconds)-t0
units = ' seconds'
if runtime gt 60 then begin
   runtime /= 60
   units = ' minutes'
   if runtime gt 60 then begin
      runtime /= 60
      units = ' hours'
      if runtime gt 24 then begin
         runtime /= 24
         units = ' days'
      endif
   endif
endif

runtime = strtrim(string(runtime,format='(f100.2)'),2)

format = '("EXOFAST_DEMC: done in ",a,a,"; took ",a,"% of trial steps")'
print, format=format, runtime, units, acceptancerate

end

    
