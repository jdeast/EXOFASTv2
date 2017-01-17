;+
; NAME:
;   EXOFASTv2
;
; PURPOSE: 
;   Simultaneously fits RV and/or transit data for multiple
;   planets, transits, RV sources. Please cite Eastman et al., 2013
;   (http://adsabs.harvard.edu/abs/2013PASP..125...83E) if you make
;   use of this routine in your research. Please report errors or bugs
;   to jeastman@lcogt.net
;
; CALLING SEQUENCE:
;   exofast, [RVPATH=, TRANPATH=, BAND=, PRIORS=, PREFIX=, /CIRCULAR,
;             /NOSLOPE, /SECONDARY, /UPDATE, PNAME=, SIGCLIP=, NTHIN=,
;             MAXSTEPS=, MINPERIOD=, MAXPERIOD=, NMIN=, /DISPLAY,
;             /DEBUG, RANDOMFUNC=, SEED=, /SPECPRIORS, /BESTONLY,
;             NINTERP=, EXPTIME=, /LONGCADENCE]
;
; INPUTS:
;
;   RVPATH      - The path to the RV data file. The file must have 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) RV (m/s)
;                   3) err (m/s)
;                 NOTE: The units must be as specified, or the fit
;                 will be wrong or fail.
;                 NOTE 2: If omitted, just the transit data will be fit
;   TRANPATH    - The path to the transit data file. The file must
;                 have at least 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) Normalized flux
;                   3) err
;                   4) Detrend parameter 1
;                   ....
;                   N+3) Detrend parameter N
;                 NOTE: The units must be as specified, or the fit
;                 will be wrong or fail.
;                 NOTE 2: If omitted, just the RV data will bit fit
;   BAND        - The bandpass of the observed transit (see quadld.pro
;                 for allowed values).
;                 NOTE: only required if TRANPATH is specified.
;   PRIORS      - Priors on each of the parameters. See EXOFAST_CHI2
;                 for parameter definitions. Must be an N x 2
;                 elements array, where N is 15 + the number of
;                 detrending variables. If set, all non-zero values of
;                 priors[0,*] are used to start the fit, and a penalty
;                 equal to:
;                    total(((pars-priors[0,*])/priors[1,*])^2)
;                 is added to the total chi^2. 
;
;                 NOTE 1: For no prior, set the width to infinity. i.e., 
;                    priors[1,*] = !values.d_infinity.
;                 NOTE 2: Typically, the default starting guesses are
;                 good enough and need not be changed.
;
;                 TRANSIT+RV or TRANSIT-ONLY FITS: Priors on Teff and
;                 [Fe/H] (priors[*,11:12]) must be specified, either by
;                 declaring the priors array or specifying a planet
;                 name and setting the /SPECPRIORS keyword.
;
;                 RV-ONLY FITS: Priors on logg, Teff and [Fe/H]
;                 (priors[*,10:12]) must be specified, either by
;                 declaring the priors array or specifying a planet
;                 name and setting the /SPECPRIORS keyword.
;
; OPTIONAL INPUTS:
;   PREFIX      - Each of the output files will have this string as a
;                 prefix. Default is RVFILE without the
;                 extension. Cannot contain an underscore, "_".
;   MINPERIOD   - The minimum period to consider. The default is 1 day.
;   MAXPERIOD   - The maximum period to consider. The default is the
;                 range of the RV input times.
;   NMIN        - The number of minima in the Lomb-Scargle Periodogram
;                 to fit a full Keplerian orbit. Default is 5. If the
;                 eccentricity is large, this may need to be
;                 increased. The execution time of the initial global
;                 fit is directly proportional to this value (though
;                 is a small fraction of the total time of the MCMC fit).
;   PNAME       - If set, starting values from exoplanets.org for the
;                 planet PNAME will be used to seed the fit. It is
;                 insensitive to extra spaces and capitalization.
;                 e.g., "WASP-12 b" is the same as "wasp-12b".
;   SIGCLIP     - If set, an iterative fit will be performed excluding
;                 data points more than SIGCLIP*error from the best
;                 fit. Default is no clipping.
;   NTHIN       - If set, only every NTHINth element will be
;                 kept. This typically doesn't affect the
;                 resultant fit because there is a high correlation
;                 between adjacent steps and has the advantage of
;                 improved memory management and faster generation of
;                 the final plots.
;   MAXSTEPS    - The maximum number of steps to take in the MCMC
;                 chain. Note that a 32-bit installation of IDL
;                 cannot allocate more than 260 million
;                 double-precision numbers, and redundant copies of
;                 each parameter are required. A very large number
;                 will cause memory management problems. Default is
;                 100,000.
;   RANDOMFUNC  - A string specifying the name of the random number
;                 generator to use. This generator must be able to
;                 return 1,2 or 3 dimensional uniform or normal random
;                 deviates. Default is 'EXOFAST_RANDOM' (which is slow
;                 but robust).
;   SEED        - The seed to the random number generator used for the
;                 MCMC fits. Be sure not to mix seeds from different
;                 generators. The default is -systime(/seconds).
;   NINTERP     - If set, the each model data point will be an average
;                 of NINTERP samples over a duration of EXPTIME
;                 centered on the input time.
;   EXPTIME     - If set, the each model data point will be an average
;                 of NINTERP samples over a duration of EXPTIME
;                 centered on the input time.
;                 For a dicussion on binning, see:
;                 http://adsabs.harvard.edu/abs/2010MNRAS.408.1758K
;   RSTAR       - A two-element array specifying the prior (rstar[0])
;                 and prior width(rstar[1]) for the stellar radius, in
;                 units of Rsun. 
;                 This will override the use of the Torres relation. 
;                 Teff and logg priors should still be specified to
;                 constrain the limb darkening.
;                 This will typically be, but need not be used
;                 in conjunction with MSTAR. This will typically be
;                 used for MSTAR < 0.6 Msun, when Torres is no longer
;                 applicable.
;   MSTAR       - Same as RSTAR, but for the stellar mass, in units of
;                 Msun.
;   MAXGR       - The maximum Gelman Rubin statistic that is
;                 considered well-mixed (default=1.01)
;   MINTZ       - The minimum number of independent draws that is
;                 considered well-mixed (default=1000)
;
; OPTIONAL KEYWORDS:
;   CIRCULAR  - If set, the fit will be forced to be circular (e=0,
;               omega_star=pi/2)
;   NOSLOPE   - If set, it will not fit a slope to the RV data.
;   SECONDARY - If set, fit a secondary eclipse. This feature is not
;               well-tested -- use at your own risk.
;   UPDATE    - Update the local copy of the exoplanets.org file (only
;               applied if PNAME is specified).
;   DISPLAY   - If set, the plots (below) will be displayed via
;               ghostview (gv).
;   SPECPRIORS- If set (and PNAME is specified), spectroscopic priors
;               from exoplanets.org will be used for logg, Teff, and
;               [Fe/H].
;   BESTONLY  - If set, only the best fit (using AMOEBA) will be
;               performed.
;   DEBUG     - If set, various debugging outputs will be printed and
;               plotted.
;   LONGCADENCE - If set, EXPTIME=29.425 and NINTERP=10 are set to handle
;                 long cadence data from Kepler.
;   YY        - If set, use the YY evolutionary tracks (see
;               http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2001ApJS..136..417Y)
;               to constrain the mass/radius relation of the star
;               instead of the Torres relation. This requires an
;               additional installation step (see
;               $EXOFAST_PATH/other/README), and significantly slows
;               down the code (~10x), but is more accurate and more
;               general than the Torres relation, and determines the
;               stellar age.
;   TIDES     - If set, when (1-Rstar/a-rp/a) < e < (1-3*Rstar/a), we
;               set the eccentricity to zero (see Eastman, 2013).
;
; OUTPUTS:
;
;   Each of the output files will be preceeded by PREFIX (defined
;   above). The first "?" will be either "c" for circular if /circular
;   is set or "e" for eccentric. The second "?" will be either "f" for
;   flat (if /noslope is set) or "m "for slope.
;
;   mcmc.?.?.idl   - An IDL save file that contains the full chains
;                    for each parameter, including derived paramters,
;                    the corresponding names, the chi2 at each link,
;                    and the index of the burn-in period.
;   pdfs.?.?.ps    - A postscript plot of each posterior distribution
;                    function, 8 to a page.
;   covars.?.?.ps  - A postscript plot of the covariances between each
;                    parameter, 16 to a page. The title of each plot is the
;                    correlation coefficient.
;   median.?.?.tex - The LaTeX source code for a deluxe table of the median
;                    values and 68% confidence interval, rounded appropriately. 
;   best.?.?.tex   - The LaTeX source code for a deluxe table of the best fit
;                    values and 68% confidence interval, rounded appropriately.
;   model.?.?.ps   - An postscript plot of the best-fit models and residuals.
;   model.?.?.rv   - A text file containing the best-fit model RV for each
;                    input time
;   model.?.?.flux - A text file containing the best-fit model flux for each
;                    input time
;
;   NOTE: To extract a single page out of a multi-page PS file, use:
;           psselect -p# input.ps output.ps
;         where # is the page number to extract.
;
; COMMON BLOCKS:
;   CHI2_BLOCK:
;     RV      - A structure with three tags, corresponding to each
;               column of RVFILE.
;     TRANSIT - A structure with 3 + N tags, corresponding to each
;               column of the TRANFILE.
;     PRIORS  - The priors on each of the parameters
;     BAND    - The band of the observed transit
;     DEBUG   - A flag to specify if the data/fit should be
;               plotted. This is intended for debugging only; it
;               dramatically increases the execution time.
;     NORV    - A boolean to indicate the RV should not be fit
;     NOTRAN  - A boolean to indicate the transit should not be fit
;   RV_BLOCK:
;     RVDATA  - A structure with three tags, corresponding to each
;               column of RVFILE. This is to prefit the data with a
;               lomb-scargle periodogram
; EXAMPLES:  
;   ;; fit HAT-P-3b example data (using priors from exoplanets.org):
;   exofast, rvpath='hat3.rv',tranpath='hat3.flux',pname='HAT-P-3b',$
;            band='Sloani',/circular,/noslope,/specpriors,minp=2.85,maxp=2.95
;
;   ;; fit HAT-P-3b without using exoplanets.org 
;   ;; include a prior on the period too 
;   ;; this reproduces results from Eastman et al., 2013
;   per = 2.899703d0
;   uper = 5.4d-5
;   priors = dblarr(2,19)
;   priors[1,*] = !values.d_infinity
;   priors[*, 3] = [alog10(per),uper/(alog(10d0)*per)] ;; logp prior
;   priors[*,10] = [4.61d0,0.05d0] ;; logg prior
;   priors[*,11] = [5185d0,  80d0] ;; Teff prior
;   priors[*,12] = [0.27d0,0.08d0] ;; logg prior
;   exofast,rvpath='hat3.rv',tranpath='hat3.flux',band='Sloani',/circular,$
;           /noslope,priors=priors,minp=2.85,maxp=2.95,prefix='example.'
;
; MODIFICATION HISTORY
; 
;  2013/05 -- Complete rewrite of exofast.pro. More general (fits
;             multiple planets, mulitple bands, multiple instrumental
;             offsets). Now easily extensible.
;-
pro exofastv2, rvpath=rvpath, tranpath=tranpath, band=band, priors=priors, $
               prefix=prefix,$
               circular=circular,fitslope=fitslope, secondary=secondary, $
               rossiter=rossiter,$
               update=update, pname=pname, sigclip=sigclip,$
               nthin=nthin, maxsteps=maxsteps, $
               minperiod=minperiod, maxperiod=maxperiod, nmin=nmin, $
               display=display,debug=debug, randomfunc=randomfunc,seed=seed,$
               specpriors=specpriors, bestonly=bestonly, plotonly=plotonly,$
               longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
               logfile=logfile, $
               maxgr=maxgr, mintz=mintz, yy=yy, tides=tides, nplanets=nplanets, $
               priorfile=priorfile, fitrv=fitrv, fittran=fittran,ttv=ttv, earth=earth

;; this is the stellar system structure
COMMON chi2_block, ss

;; name of the chi square function
chi2func = 'exofast_chi2v2'

;; compile all routines now to keep output legible 
;; resolve_all doesn't interpret execute; it's also broken prior to IDL v6.4(?)
if double(!version.release) ge 6.4d0 then $
   resolve_all, resolve_function=[chi2func,'exofast_random'],/cont,/quiet

;; create the master structure
ss = mkss(rvpath=rvpath, tranpath=tranpath, nplanets=nplanets, $
          debug=debug, priorfile=priorfile, fitrv=fitrv, fittran=fittran, $
          circular=circular,fitslope=fitslope, fitquad=fitquad,ttv=ttv, $
          rossiter=rossiter,longcadence=longcadence, earth=earth)

;mcmc2str, pars, ss
;stop

;; derive the step parameters from the priors
ok = pars2step(ss)

;; output to file instead of screen
if n_elements(logfile) ne 0 then openw, lun, logfile, /get_lun $
else lun = -1
if n_elements(mintz) eq 0 then mintz = 1000d0
if n_elements(maxgr) eq 0 then maxgr = 1.01d0
if n_elements(maxsteps) eq 0 then maxsteps = 100000L
if n_elements(nthin) eq 0 then nthin = 1L
if n_elements(sigclip) ne 1 then sigclip = !values.d_infinity ;; no clipping
if n_elements(nmin) eq 0 then nmin=5
;; use robust (slower) random number generator by default
if n_elements(randomfunc) eq 0 then randomfunc = 'exofast_random'

;; default prefix for all output files (filename without extension)
if n_elements(prefix) eq 0 then prefix = 'planet.'
basename = file_basename(prefix)

;; Constants

;; 2014 CODATA 
;; http://arxiv.org/abs/1507.07956
GSI = 6.67408d-11      ;; m^3/(kg*s^2) 


;; IAU Resolution B3, 2015 
;; http://arxiv.org/abs/1510.07674
rsun = 6.957d8         ;; meters
GMsun = 1.3271244d20   ;; m^3/s^2
day = 86400d0          ;; seconds
G = GMsun/rsun^3*day^2 ;; R_sun^3/(m_sun*day^2)
msun = GMsun/GSI       ;; kg

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010

;; plot the data + starting guess
if keyword_set(plotonly) then begin
   modelfile = prefix + 'model.ps'
   bestchi2 = call_function(chi2func, psname=modelfile, $
                            modelrv=modelrv, modelflux=modelflux)
   if keyword_set(display) then spawn, 'gv ' + modelfile + ' &'
   if lun ne -1 then free_lun, lun
   return
endif

pars = str2pars(ss,scale=scale,name=name)
;; these are the starting values for all step parameters
forprint, indgen(n_elements(name)), name, pars, scale,/t,format='(i3,x,a15,x,f14.6,x,f14.6)'

best = exofast_amoeba(1d-8,function_name=chi2func,p0=pars,scale=scale,nmax=1d5)
if best[0] eq -1 then begin
   printf, lun, 'ERROR: Could not find best combined fit'
   if lun ne -1 then free_lun, lun
   return
endif

;; output the best-fit model fluxes/rvs
bestchi2 = call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux)
plottran, ss, psname=prefix + 'model.ps'

;; do the MCMC fit
if not keyword_set(bestonly) then begin
;ss.debug=1
   exofast_demc, best, chi2func, pars, chi2=chi2,$
                 nthin=nthin,maxsteps=maxsteps,$
                 burnndx=burnndx, seed=seed, randomfunc=randomfunc0, $
                 gelmanrubin=gelmanrubin, tz=tz, maxgr=maxgr, mintz=mintz, $
                 derived=rstar
   if pars[0] eq -1 then begin
      printf, lun, 'MCMC Failed to find a stepping scale. This usually means one or more parameters are unconstrained by the data or priors.'
      if lun ne -1 then free_lun, lun
      return
   endif 

   bad = where(tz lt mintz or gelmanrubin gt maxgr,nbad)
   if bad[0] ne -1 then begin
      printf, lun, 'WARNING: The Gelman-Rubin statistic indicates ' + $
              'the following parameters are not well-mixed'
      printf, lun, '    Parameter   Rz     Tz'
      for i=0, nbad-1 do printf, lun, name[bad[i]], gelmanrubin[bad[i]],tz[bad[i]], format='(a13,x,2(f0.4,x))'
   endif
   print, 'Synthesizing results; for long chains and/or many fitted parameters, this may take up to 15 minutes'

   ;; combine all chains
   sz = size(pars)
   npars = sz[1]
   nsteps = sz[2]
   nchains = sz[3]
   pars = reform(pars,npars,nsteps*nchains)
   chi2 = reform(chi2,nsteps*nchains)
   minchi2 = min(chi2,bestndx)

   rstar = reform(rstar,1,nsteps*nchains)
endif else begin
   pars = reform(best[tofit],n_elements(tofit),1)
   bestndx = 0
endelse

;; generate the model fit from the best MCMC values, not AMOEBA
bestamoeba = best
best = pars[*,bestndx]
modelfile = prefix + 'model.ps'
;bestchi2 = call_function(chi2func,best,psname=modelfile, $
;                         modelrv=modelrv, modelflux=modelflux)

;; make a new stellar system structure with only fitted and derived
;; parameters, populated by the pars array
;mcmcss = mcmc2str(pars, ss)


;; derive all parameters where derive=1 (final values)
mcmcss = mkss(rvpath=rvpath, tranpath=tranpath, nplanets=nplanets, $
              debug=debug, priorfile=priorfile, fitrv=fitrv, fittran=fittran, $
              circular=circular,fitslope=fitslope, fitquad=fitquad,$
              nvalues=nsteps*nchains,ttvs=ttv,longcadence=longcadence)
mcmcss.burnndx = burnndx
mcmcss.star.rstar.value = rstar

pars2str, pars, mcmcss

;; derive all parameters
derivepars, mcmcss

;; save the chains for additional analysis
idlfile = prefix + 'mcmc.idl'
save, mcmcss, filename=idlfile

;; output filenames
label = "tab:" + basename
caption = "Median values and 68\% confidence interval for " + basename
parfile = prefix + 'pdf.ps'
covarfile = prefix + 'covar.ps'
texfile = prefix + 'median.tex'

exofast_plotdist2, mcmcss, pdfname=parfile, covarname=covarfile;,/nocovar
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

;; display all the plots, if desired
if keyword_set(display) then begin
   spawn, 'gv ' + parfile + ' &'
   spawn, 'gv ' + covarfile + ' &'
   spawn, 'gv ' + modelfile + ' &'
endif

if lun ne -1 then free_lun, lun

end
