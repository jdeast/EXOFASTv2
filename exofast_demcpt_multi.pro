pro exofast_callback, olddet, oldpars, oldchi2, oldderived, nderived, nswap, naccept, temps, randomfunc, swapped, $
                      det=det, dpar=dpar, newchi2=newchi2, fac=fac, j=j, m=m, k=k, newpars=newpars, thread=thread

  if n_elements(thread) ne 0 then begin
     det = thread.obridge->GetVar('det')
     if nderived gt 0 then dpar = thread.obridge->GetVar('dpar')
     newchi2 = thread.obridge->GetVar('newchi2')
     newpars = thread.newpars
     fac = thread.fac
     j = thread.j
     k = thread.k
     m = thread.m
     thread.status=0
     thread.j = -1L
     thread.k = -1L
     thread.m = -1L

     ;print, 'thread ' + strtrim(j,2) + ' ' + strtrim(m,2) + ' took ' + strtrim(systime(/seconds)-thread.start,2)

  endif

  C = fac*olddet[j,m]/det*exp(temps[m]*(oldchi2[j,m] - newchi2)/2d0)

  if m eq 0 then begin
;     print, j, m, c, newchi2, oldchi2[j,m]
;     print, oldpars[*,j,m]
;     print, newpars
  endif

  ;; accept the step; update values
  if call_function(randomfunc,seed) lt C then begin
     if swapped then nswap++ $
     else naccept++
     olddet[j,m] = det
     oldpars[*,j,m] = newpars
     oldchi2[j,m] = newchi2
     if nderived gt 0 then oldderived[*,j,m] = dpar

     if m eq 0 then begin
;        if swapped then print, 'swapped!' $
;        else stop;print, 'accepted!'
     endif                           

  endif ;; else keep previous values
  
;  if m eq 0 then print

;  if m eq 0 and not swapped then stop

;stop

end
  
;+
; NAME:
;   exofast_demc
; PURPOSE:
;   Make NCHAINS independent Markov Chain Monte Carlo chains to fit data
;   and return the parameter distributions.
;
; DESCRIPTION:
;   Begins by determining the correct stepping scale by finding
;   differences in each parameters that yield a delta chi^2 = 1, then
;   uses 5*scale*randomn offsets in each parameter for the starting
;   points of each chain.
;
;   It then begins a Differential Evolution Markov Chain Monte Carlo
;   fit (ter Braak, 2006)
;   http://www.stat.columbia.edu/~gelman/stuff_for_blog/cajo.pdf.  The
;   only slight modification to this basic alogorithm is that we
;   determine the magnitude of the uniform random deviate for each
;   parameter. We have found the dynamic range of our parameters can
;   be large, and thus not well suited to the one-size fits all
;   approach recommended there.
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
;   useable step so as not to be biased by the starting criteria.
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
; SYSTEM VARIABLES:
;   !STOPNOW   - If set, will break the MCMC loop at the next complete
;                step. It must finish the current step, including all
;                NTHIN and NCHAIN evaluations. To set it, type
;                <CONTROL> + C
;                !STOPNOW = 1 
;                .con
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
;-

pro exofast_demcpt_multi, bestpars,chi2func,pars,chi2=chi2, tofit=tofit,$
                    scale=scale,seed=seed,randomfunc=randomfunc,$
                    nthin=nthin, maxsteps=maxsteps, maxtime=maxtime, ntemps=ntemps, tf=tf, $
                    keephot=keephot, hotpars=hotpars, allhotpars=allhotpars, hotchi2=hotchi2, $
                    dontstop=dontstop,$
                    nchains=nchains, angular=angular,$
                    burnndx=burnndx, removeburn=removeburn, goodchains=goodchains,$
                    gelmanrubin=gelmanrubin, tz=tz, maxgr=maxgr, mintz=mintz, $
                    derived=derived, resumendx=resumendx,stretch=stretch, logname=logname, acceptancerate=acceptancerate,$
                    thread_array=thread_array

COMMON chi2_block, ss

;; default values
if n_elements(maxsteps) eq 0 then maxsteps = 1d5
if n_elements(maxtime) eq 0 then maxtime = !values.d_infinity
if n_elements(nthin) eq 0 then nthin = 1L
if n_elements(tofit) eq 0 then tofit = indgen(n_elements(bestpars))
if n_elements(mintz) eq 0 then mintz = 1000d0
if n_elements(maxgr) eq 0 then maxgr = 1.01d0
if n_elements(loglun) eq 0 then loglun = -1
if n_elements(ntemps) eq 0 then ntemps=1
lasttz = 0
lastgr = 10

;; the number of parallel temperatures to run
if ntemps eq 1 then begin
   temps = 1d0
   ntemps = 1d0   
   swapped = 0B
endif else begin
   if n_elements(tf) eq 0 then tf = 200d0
   temps = (1d0/tf)^(dindgen(ntemps)/(ntemps-1d0))
endelse

defsysv, '!STOPNOW', 0B
nfit = n_elements(tofit)
!except = 0 ;; don't display errors for NaN and infinity

;; defaults for Differential Evolution MCMC
gamma = 2.38d0/sqrt(2d0*nfit)
if n_elements(nchains) eq 0 then nchains=(2d0*nfit > 3)

;; prior to IDL v8.2, randomu used ran1, which is terrible
;; later versions use Mersenne Twister, which is ok for MCMC
;; exofast_random uses ran3, which is even better
;;    but the native IDL implmentation is slow
if n_elements(randomfunc) eq 0 then begin
   if !version.release ge 8.2 then randomfunc = 'randomu' $
   else randomfunc = 'exofast_random'
endif


;; angular must be indexed within tofit to use with gelman-rubin
nang = n_elements(angular)
if nang gt 0 then gelmanangular=intarr(nang)
for i=0, nang-1 do gelmanangular[i] = (where(tofit eq angular[i]))(0)

junk = call_function(chi2func, bestpars, determinant=det, derived=dpar)
nderived = n_elements(dpar)

oldpars = dblarr(nfit,nchains,ntemps) ;; parameters in the previous step
oldchi2 = dblarr(nchains,ntemps)       ;; chi2 of the previous step
olddet = dblarr(nchains,ntemps)        ;; determinant of the previous step
if keyword_set(keephot) then begin
   hotpars = dblarr(nfit,maxsteps,nchains,ntemps-1)
   hotchi2 = dblarr(maxsteps,nchains,ntemps-1)
endif

;; initialize arrays
if n_elements(resumendx) ne 0 then begin
   ;; this is a(n untested) placeholder for a future addition to be
   ;; able to restart runs from saved IDL files

   sz = size(pars)
   nfit = sz[1]
   maxsteps = sz[2]
   nchains = sz[3]
   newpars = pars[*,resumendx,0]
   oldpars[*,*,0] = pars[*,resumendx-1,0]
   if nderived gt 0 and n_elements(derived) ne 0 then begin
      printandlog, "ERROR: there are derived parameters; you must supply a DERIVED array to resume from previous run", logname
      stop
   endif
   
   if ntemps gt 0 then oldpars[*,*,1:n_elements(oldpars[0,0,*])-1] = hotpars
   
   for j=0, nchains-1 do begin
      for m=0, ntemps-1 do begin
         junk = call_function(chi2func, pars[*,resumendx-1,j], determinant=det, derived=dpar)
         olddet[j,m] = det
      endfor
   endfor

endif else begin
   ;; this is the expected/tested switch here
   resumendx = 1L
   pars = dblarr(nfit,maxsteps,nchains)
   chi2 = dblarr(maxsteps,nchains)
   newpars = bestpars
   oldpars[*,0,*] = replicate(1d0,ntemps)##bestpars

   ;; get the MCMC step scale for each parameter
   if n_elements(scale) eq 0 then $
      scale = exofast_getmcmcscale(bestpars,chi2func,tofit=tofit,angular=angular,/skipiter,logname=logname)
   if scale[0] eq -1 then begin
      printandlog, 'No scale found',logname
      pars = -1
      return
   endif

   printandlog, 'Beginning chain initialization', logname
   ;; initialize each chain
   for j=0, nchains-1 do begin
      for m=0, ntemps-1 do begin
         ;; repeat until a finite chi^2 
         ;; i.e., don't start outside a boundary; you'll never get back
         niter = 0d0
         repeat begin
            ;; start 5 steps from best value. If that's not allowed
            ;; (infinite chi^2), slowly approach 0 steps away from the best value
            ;; niter should never be larger than ~5000
;            if j eq 0 then factor = 0d0 else factor = (sqrt(500d0/nfit) > 3d0) ;; first chain starts at best fit
            if j eq 0 then factor = 0d0 else factor = (sqrt(500d0/nfit) < 3d0) ;; first chain starts at best fit
            oldpars[*,j,m] = bestpars[tofit] + $
                             factor/exp(niter/1000d)*scale*call_function(randomfunc,seed,nfit,/normal)
            
            newpars[tofit] = oldpars[*,j,m]
            ;; find the chi^2
            oldchi2[j,m] = call_function(chi2func, newpars, determinant=det, derived=dpar)
            niter += 1
         endrep until finite(oldchi2[j,m])

         if n_elements(dpar) ne 0 then oldderived[*,j,m] = dpar
         olddet[j,m] = det
         oldpars[*,j,m] = newpars[tofit]
      endfor
      chi2[0,j] = oldchi2[j,0] ;; populate the first element of the chi^2 array
   endfor

   printandlog, 'Done initializing chains', logname
   printandlog, '', logname

   if n_elements(dpar) ne 0 then oldderived = reform(derived[*,0,*],nderived,nchains)

endelse

pars[0:nfit-1,resumendx-1,*] = oldpars[*,*,0]
if n_elements(hotpars) gt 0 then begin
   hotpars[*,0,*,*] = oldpars[*,*,1:ntemps-1]
   hotchi2[0,*,*] = oldchi2[*,1:ntemps-1]
endif      


a = 2d0 ;; for affine invariant step
nextrecalc = 100L > (resumendx/0.9d0) ;; when to recalculate the convergence stats
npass = 1L ;; the number of consecutive times it has passed the convergence criteria
nstop = 0L ;; the step it was stopped at (either converged or interrupted)
naccept = 1d0 ;; the number of accepted steps
tz0 = 0d0 ;; 
tzsteps = 0L
alreadywarned = 0L
t0 = systime(/seconds)
nswap = 0L
nswapattempt = 0L
nattempt = 0L
defsysv, '!GDL', exists=runninggdl
t00 = systime(/seconds)
nthreads=n_elements(thread_array)

;; start MCMC chain
for i=resumendx,maxsteps-1L do begin

;   if i ne resumendx then printandlog, 'the ' + strtrim(i,2) + 'th step took ' + strtrim(systime(/seconds)-t00,2) + ' seconds', logname
;   t00 = systime(/seconds)

   ;; automatically thin the chain (saves memory)
   for k=0L, nthin-1L do begin

      threadndx = -1L
      for j=0L, nchains-1L do begin
      
         for m=0L, ntemps-1L do begin

            ;; Parallel Tempering step (no need to recompute the chi2)
            ;; http://www.lindonslog.com/programming/openmp/parallel-tempering-algorithm-c/
            ;; attempt swap between next hottest chain with probability 1/2
            if m ne ntemps-1 and call_function(randomfunc,seed) lt 0.5d0 then begin
               newpars[tofit] = oldpars[*,j,m+1]
               newchi2 = oldchi2[j,m+1]
               det = olddet[j,m+1]
               if nderived gt 0 then dpar = oldderived[*,j,m+1]             
               nswapattempt++
               swapped = 1B
               fac = 1d0
               exofast_callback, olddet, oldpars, oldchi2, oldderived, nderived, $
                                 nswap, naccept, temps, randomfunc, swapped, $
                                 det=det, dpar=dpar, newchi2=newchi2, fac=fac, $
                                 j=j, m=m, newpars=newpars
            endif else begin

               ;; a random chain
               r1=floor(call_function(randomfunc,seed)*(nchains-1))
               if r1 eq j then r1 = long(nchains)-1L ;; can't use self, reassign to last
               
               if keyword_set(stretch) then begin
                  ;; affine invariant "stretch move" step 
                  ;; see Goodman & Weare, 2010 (G10), eq 7 
                  ;; http://msp.org/camcos/2010/5-1/camcos-v5-n1-p04-p.pdf
                  
                  ;; a random number between 1/a and a, weighted toward lower
                  ;; values to satisfy the symmetry condition (G10, eq 8)
                  z = ((a-1d0)*call_function(randomfunc,seed) + 1d0)^2d0/2d0 
                  
                  ;; use the most recent step (Goodman, pg 69)
                  ;; define the new step (eq 7 in G10)
                  newpars[tofit] = oldpars[*,r1,m] + z*(oldpars[*,j,m]-oldpars[*,r1,m])
                  fac = z^(nfit-1) ;; GW10, last equation, pg 70
               endif else begin
                  ;; differential evolution mcmc step (Ter Braak, 2006)
                  
                  ;; pick a second random chain (that isn't self or the first)
                  r2=floor(call_function(randomfunc,seed)*(nchains-2))
                  if r2 ge r1 then r2+=1
                  if r2 ge j then r2+=1
                  if r2 eq r1 then r2+=1
                  
                  newpars[tofit] = oldpars[*,j,m] + $ 
                                   gamma*(oldpars[*,r1,m]-oldpars[*,r2,m] + $
                                          (call_function(randomfunc,seed,nfit)-0.5d0)*scale/10d0)
                  
                  fac = 1d0

               endelse 
               
               ;; calculate the chi^2 of the new step
               ;; this is super expensive in general
               ;; farm this out to multiple threads 
               if nthreads gt 1 then begin
                  ;; find a thread that's not executing
                  breakloop = 0B
                  repeat begin
                     for ithread=threadndx+1,nthreads+threadndx do begin   
                        thread = ithread mod nthreads ;; check the oldest thread first

                        ;; Status values: 0=idle, 1=Executing Command, 2=Completed, 3=Error, 4=Aborted 
                        thread_status = thread_array[thread].obridge->status()

                        ;print, thread, thread_status

                        ;; if this thread is not executing, use it
                        if thread_status ne 1 then begin 
                           
                           ;; this thread is not executing, use it
                           ;; (and clean up if necessary)
                           breakloop = 1B ;; break repeat loop
                           threadndx = thread
                           
                           ;; clean up the thread if it's done, not just idle
                           if thread_array[thread].status ne 0 then begin
                              swapped = 0B

                              thisthread = thread_array[thread]
                              exofast_callback, olddet, oldpars, oldchi2, oldderived, $
                                                nderived, nswap, naccept, temps, randomfunc, $
                                                swapped, thread=thisthread
                              thread_array[thread] = thisthread

                           endif                           
                           break ;; end for loop                           
                        endif
                     endfor
                  endrep until breakloop

                  ;; now execute the new thread
                  thread_array[threadndx].start = systime(/seconds)
                  thread_array[threadndx].obridge->setvar, 'newpars', newpars
                  thread_array[threadndx].obridge->execute, 'newchi2 = call_function(chi2func,newpars,det=det,derived=dpar)',/nowait
                  thread_array[threadndx].k = k
                  thread_array[threadndx].j = j
                  thread_array[threadndx].m = m
                  thread_array[threadndx].fac = fac
                  thread_array[threadndx].newpars = newpars
                  thread_array[threadndx].status = 1

;                  thread_status = thread_array[threadndx].obridge->status()
;                  newchi2 = thread_array[threadndx].obridge->GetVar('newchi2')
;                  newchi22 = call_function(chi2func,newpars,det=det,derived=dpar)
;                  exofast_callback, olddet, oldpars, oldchi2, oldderived, nderived, $
;                                    nswap, naccept, temps, randomfunc, swapped, det=det, $
;                                    dpar=dpar, newchi2=newchi2, fac=fac, j=j, m=m, newpars=newpars


;                  print, newchi22, newchi2, thread_status
               endif else begin

                  ;; if nthreads=1, do it serially
                  swapped = 0B
                  newchi2 = call_function(chi2func,newpars,det=det,derived=dpar)
                  exofast_callback, olddet, oldpars, oldchi2, oldderived, nderived, $
                                    nswap, naccept, temps, randomfunc, swapped, det=det, $
                                    dpar=dpar, newchi2=newchi2, fac=fac, j=j, m=m, newpars=newpars

;                  print, newpars
;                  print, oldpars[*,j]
;                  print, newchi2, oldchi2[j]
;stop

               endelse
               nattempt++
            endelse
            
         endfor ;; each temperature

         ;; update the link in the chain with the temp=1 value
         pars[*,i,j] = oldpars[*,j,0]
         chi2[i,j] = oldchi2[j,0]
         if n_elements(dpar) ne 0 then derived[*,i,j] = oldderived[*,j,0]
         
         if n_elements(hotpars) gt 0 then begin
            hotpars[*,i,j,*] = oldpars[*,j,1:ntemps-1]
            hotchi2[i,j,*] = oldchi2[j,1:ntemps-1]
         endif      

      endfor ;; each chain

      ;; wait for all threads to finish the j chains and m temps 
      ;; of the kth thin and ith step
      if nthreads gt 1 then begin
         repeat begin
            breakloop = 1B
            for thread=0L,nthreads-1 do begin
               status = thread_array[thread].obridge->status()
               ;; Status values:  0=idle, 1=Executing Command, 2=Completed, 3=Error,4 = Aborted 
               if status ne 1 then begin
                  if thread_array[thread].status eq 1 then begin
                     ;; thread is ready for cleanup
                     thisthread = thread_array[thread]
                     swapped = 0B ;; swapped steps are cheap (not threaded)
                     exofast_callback, olddet, oldpars, oldchi2, oldderived, nderived, $
                                       nswap, naccept, temps, randomfunc, swapped, thread=thisthread
                     thread_array[thread] = thisthread
                  endif
               endif else breakloop = 0B ;; at least one thread is still active, keep waiting
            endfor
         endrep until breakloop
      endif
      
   endfor ;; nthin
   
   ;; stop the fit if it has exceeded the maximum allowed runtime
   tnow = systime(/seconds)
   if tnow-t0 gt maxtime then begin
      printandlog, '', logname
      printandlog, 'Time limit reached, stopping the run', logname
      !stopnow=1
   endif

   ;; this allows the user to break at any point
   defsysv, '!STOPNOW', exists=stopnowdefined
   if stopnowdefined then begin
      if !STOPNOW EQ 1 then break
   endif

   ;; Test for convergence as outlined in Ford 2006
   ;; must be converged for 6 consecutive passes
   ;; tz > 1000 and Gelman-Rubin < 1.01 => converged 
   if i eq nextrecalc then begin
      ;; discard the burn-in
      burnndx = getburnndx(chi2[0:i,*],goodchains=goodchains,/silent)

      ;; calculate the Gelman-Rubin statistic (remove burn-in)
      converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:i,goodchains],$
                                      gelmanrubin,tz,angular=gelmanangular,$
                                      mintz=mintz,maxgr=maxgr)
      lasttz = min(tz)
      lastgr = max(gelmanrubin)

      ;; estimate the number of steps it will take until it is converged
      tz0 = [tz0,min(tz)]
      tzsteps = [tzsteps,i]
      ntz = n_elements(tz0)
      
      ;; tz starts around the number of chains, even though the mixing
      ;; is poor. Don't start counting until it's double that.
      good = where(tz0 gt nchains*2d0,ntz)
      if ntz gt 0 then ntz = n_elements(tzsteps[good[0]:*])

;      ;; some debugging plots
;      if ntz ge 3 then begin
;         window, 0, retain=2
;         plot, tzsteps[good[0]:*], tz0[good[0]:*], psym=1,/ynoz
;         oplot, tzsteps[good[0]:*], tz0[good[0]:*]
;         coeffs = poly_fit(tzsteps[good[0]:*],tz0[good[0]:*],1,yfit=yfit)
;         oplot, tzsteps[good[0]:*], yfit, color='0000ff'x
;         stepstoconvergence = (mintz - coeffs[0])/coeffs[1]
;         print
;         print, stepstoconvergence
;      endif

      if not alreadywarned then begin
         if i gt maxsteps/20d0 and ntz gt 3 then begin
            ;; fit a line to the number of independent draws
            coeffs = poly_fit(tzsteps[good[0]:*],tz0[good[0]:*],1,yfit=fit)
            
            ;; extrapolate to 2*MINTZ, when it'll be considered converged 
            ;; MINTZ is considered converged, but this is far from exact...
            stepstoconvergence = (2d0*mintz - coeffs[0])/coeffs[1]
            
            ;; if it won't be converged, warn the user
            if stepstoconvergence gt maxsteps then begin
               bestnthin = round(stepstoconvergence*nthin/maxsteps,/l64)
               printandlog,'WARNING: The chain is not expected to be well-mixed; '+$
                       'set NTHIN to >~ ' + strtrim(bestnthin,2) + ' to fix', logname
            endif else if stepstoconvergence le 0 then begin
               printandlog, 'WARNING: The chains are not mixing. This may be a while.', logname
            endif else printandlog, string(stepstoconvergence*50d0/maxsteps,$
                                              format='("EXOFAST_DEMC: The chain is expected to be ' + $
                                              'well-mixed after ",i2,"% completed     ")'), logname
            alreadywarned = 1 ;; only display message once
         endif
      endif
      
      ;; must pass 6 consecutive times
      if converged then begin
         if nstop eq 0 then nstop = i
         nstop = i
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
   acceptancerate = strtrim(string(naccept/double(nattempt)*100,$
                                   format='(f0.2)'),2)
   swaprate = strtrim(string(nswap/double(nswapattempt)*100,$
                             format='(f6.2)'),2)

   timeleft = ((tnow-t0)*(maxsteps/(i+1d)-1d)) < (maxtime - (tnow-t0))
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
if 1 then begin
;   if !version.os_family ne 'Windows' or (i eq resumendx) then begin
;   if !version.os_family ne 'Windows' or $
;      ((i+1) mod round(maxsteps/1000d0) eq 0) or (i eq resumendx) then begin

      if ntemps gt 1 then begin
         format='("EXOFAST_DEMC: ",f0.2,"% done; acceptance rate = ",a,"%; swap rate = ",a,"%; ' + $ 
                'tz = ",f0.2," (>", a,"); GelmanRubin = ",f0.4," (<",f0.2,"); time left: ",f0.2,a,$,%"\r")'
         print, 100.d0*(i+1)/maxsteps,acceptancerate,swaprate,lasttz,strtrim(fix(mintz),2),lastgr,maxgr,timeleft,units, format=format
         
         ;; print this message to the log every 5% of the way
         if i eq resumendx or (i+1) mod round(maxsteps/20) eq 0 then begin
            format='("EXOFAST_DEMC: ",f0.2,"% done; acceptance rate = ",a,"%; swap rate = ",a,"%; ' + $ 
                   'tz = ",f0.2," (>", a,"); GelmanRubin = ",f0.4," (<",f0.2,"); time left: ",f0.2,a)'         
            printandlog, string(100.d0*(i+1)/maxsteps,acceptancerate,swaprate,lasttz,strtrim(fix(mintz),2),lastgr,maxgr,timeleft,units,format=format), logname
         endif
      endif else begin
         format='("EXOFAST_DEMC: ",f0.2,"% done; acceptance rate = ",a,"%; ' + $ 
                'tz = ",f0.2," (>", a,"); GelmanRubin = ",f0.4," (<",f0.2,"); time left: ",f0.2,a,$,%"\r")'
         print, 100.d0*(i+1)/maxsteps,acceptancerate,lasttz,strtrim(fix(mintz),2),lastgr,maxgr,timeleft,units, format=format
         
         ;; print this message to the log every 5% of the way
         if i eq resumendx or (i+1) mod round(maxsteps/20) eq 0 then begin
            format='("EXOFAST_DEMC: ",f0.2,"% done; acceptance rate = ",a,"%; ' + $ 
                   'tz = ",f0.2," (>", a,"); GelmanRubin = ",f0.4," (<",f0.2,"); time left: ",f0.2,a)'         
            printandlog, string(100.d0*(i+1)/maxsteps,acceptancerate,lasttz,strtrim(fix(mintz),2),lastgr,maxgr,timeleft,units,format=format), logname
         endif
      endelse
   endif
   
endfor
printandlog, '', logname ;; don't overwrite the final line
printandlog, '', logname ;; now just add a space

;; JDE 4/1/2022 -- this isn't right. 
;; if nstop < maxsteps (typical), it breaks the loop without incrementing i
;; this line would remove the last step in most cases
;nstop = i-1L 
nstop = (i < (maxsteps-1))

burnndx = getburnndx(chi2[0:nstop,*],goodchains=goodchains, logname=logname)
ngood = n_elements(goodchains)

;; if there are bad chains, is it better without them?
;; (it usually is unless it's a very short run)
if ngood ne nchains or npass ne 6 then begin

   ;; are all the chains mixed?
   converged = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,*],$
                                   gelmanrubin1,tz1,angular=gelmanangular, $
                                   maxgr=maxgr, mintz=mintz)
   
   ;; ngood is required to be > 2....
   if ngood gt 2 then begin
      converged2 = exofast_gelmanrubin(pars[0:nfit-1,burnndx:nstop,goodchains],$
                                       gelmanrubin2,tz2,angular=gelmanangular, $
                                       maxgr=maxgr, mintz=mintz)
      bad = where(tz2 lt mintz or gelmanrubin2 gt maxgr)
      if bad[0] ne -1 then extrabad = where(tz2[bad] lt tz1[bad] or $
                                            gelmanrubin2[bad] gt gelmanrubin1[bad]) $
      else extrabad = -1
      
   endif else begin
      ;; this should never execute
      bad = where(tz1 lt mintz or gelmanrubin1 gt maxgr)
      extrabad = -1
      printandlog, 'There were less than three bad chains, which should not be possible.',logname 
   endelse
   
   if bad[0] eq -1 then begin
      if npass ne 6 then printandlog, 'Chains are marginally well-mixed',logname $
      else begin
         printandlog,'Some chains never got below the median chi^2 and have been discarded before inference.',logname
         printandlog,'Starting a new fit at the best-fit found here may remove this message.',logname
      endelse
   endif else if extrabad[0] eq -1 and ngood ne nchains and ngood gt 2 then begin
      ;; without errant chains, it's better
      printandlog, 'Discarding ' + strtrim(round(nchains-ngood),2) + '/'+$
              strtrim(round(nchains),2) + ' chains stuck in local minima.', logname
      pars = pars[*,*,goodchains]
      chi2 = chi2[*,goodchains]
      if n_elements(derived) ne 0 then derived = derived[*,*,goodchains]
      nchains = ngood
      tz = tz2
      gelmanrubin = gelmanrubin2
      bad = where(tz lt mintz or gelmanrubin gt maxgr)
      if not converged2 then begin
         printandlog, 'Remaining chains are still not mixed. Must run longer.', logname
         printandlog, '       Parameter  Gelman-Rubin Independent Draws',logname
         for i=0L, n_elements(bad)-1 do printandlog, string(tofit[bad[i]], gelmanrubin2[bad[i]],tz2[bad[i]]),logname
      endif
   endif else begin
      printandlog, 'No obviously-errant chains. Must run longer.', logname
      printandlog, '       Parameter  Gelman-Rubin Independent Draws',logname
      for i=0L, n_elements(bad)-1 do printandlog, string(tofit[bad[i]], gelmanrubin1[bad[i]],tz1[bad[i]]), logname
   endelse
   if n_elements(gelmanrubin) eq 0 then gelmanrubin=gelmanrubin1
   if n_elements(tz) eq 0 then tz=tz1
endif

;; remove the uncalculated parameters if terminated early
;; we keep the burn in by default so people can do their own analysis
if keyword_set(removeburn) then begin
   ;; but remove the burn-in if /REMOVEBURN is set
   pars = pars[*,burnndx:nstop,*]
   chi2 = chi2[burnndx:nstop,*]
   if n_elements(derived) ne 0 then derived = derived[*,burnndx:nstop,*]
   if n_elements(hotpars) ne 0 then hotpars = hotpars[*,burnndx:nstop,*,*]
   if n_elements(hotchi2) ne 0 then hotchi2 = hotchi2[burnndx:nstop,*,*]
endif else begin
   pars = pars[*,0:nstop,*]
   chi2 = chi2[0:nstop,*]
   if n_elements(derived) ne 0 then derived = derived[*,0:nstop,*]
   if n_elements(hotpars) ne 0 then hotpars = hotpars[*,0:nstop,*,*]
   if n_elements(hotchi2) ne 0 then hotchi2 = hotchi2[0:nstop,*,*]
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

format = '("EXOFAST_DEMC: done in ",a,a,"; accepted ",a,"% of trial steps")'
printandlog, string(runtime, units, acceptancerate, format=format), logname

end

    
