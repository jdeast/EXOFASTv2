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
;   exofast, [RVPATH=, TRANPATH=, PRIORFILE=, FLUXFILE=, PREFIX=, /CIRCULAR,
;             /NOSLOPE, /SECONDARY, /UPDATE, PNAME=, SIGCLIP=, NTHIN=,
;             MAXSTEPS=, MINPERIOD=, MAXPERIOD=, NMIN=, /DISPLAY,
;             /DEBUG, RANDOMFUNC=, SEED=, /SPECPRIORS, /BESTONLY,
;             NINTERP=, EXPTIME=, /LONGCADENCE]
;
; INPUTS:
;
;   PRIORFILE  - An ASCII file with each line containing three white
;                space delimited columns: NAME, VALUE, WIDTH. NAME
;                must match a parameter.label defined in MKSS.pro. If
;                the parameter is in an array (e.g., of planets), add
;                "_i", where "i" is the zero-indexed index into the
;                array. If WIDTH is set to 0, the parameter is fixed
;                at VALUE (this is generally not recommended; it's far
;                better to apply a realistic prior). If WIDTH is set
;                to -1, it only overrides the default starting guess
;                for that parameter, but there is no penalty if the
;                model deviates from VALUE. If WIDTH is positive, a
;                gaussian prior is applied. That is, a chi^2 penalty
;                equal to ((parameter - VALUE)/WIDTH)^2 is applied to
;                the likelihood function. Here is the contents of a
;                sample priorfile for the EPXXXXXXXXX system published
;                in Eastman et al, 2017:
;
;                teff 6167 78 # Gaussian prior on T_eff of 6167 +/- 78 K
;                feh -0.04 0.1 # Gaussian prior on [Fe/H] of -0.04 +/- 0.1 dex
;                logg 4.22 0.06 # Gaussian prior on logg of 4.22 +/- 0.06
;                vsini 9400 300 # Gaussian prior on vsini of 9400 +/- 300 (ignored since vsini is not fitted without RM data)
;                # parallax 4.034 0.93 # Gaussian prior on parallax of 4.036 +/- 0.93 mas. However, it is commented out, so not applied.
;                tc_0 2457166.0543 -1 # start the fit with a TC for planet 0 (EPXXXXXXXXb) at BJD_TDB=2457166.0543
;                p_0 0.020 -1 # start the fit with Rp/Rstar for planet 0 at 0.020
;                period_0 26.847 -1 # start the fit with period for planet 0 at 26.847 days
;                tc_1 2457213.5328 -1 # start the fit with a TC for planet 1 (EPXXXXXXXXc) at BJD_TDB=2457213.5328
;                p_1 0.02 -1 # start the fit with Rp/Rstar for planet 1 at 0.020
;                period_1 39.5419 -1 # start the fit with period for planet 0 at 39.5419 days
;                tc_2 2457191.8894 -1 # start the fit with a TC for planet 2 (EPXXXXXXXXd) at BJD_TDB=2457191.8894
;                p_2 0.020 -1 # start the fit with Rp/Rstar for planet 2 at 0.020
;                period_2 125 -1 # start the fit with period for planet 2 at 125 days
;                tc_3 2457170.4990 -1 # start the fit with a TC for planet 3 (EPXXXXXXXXe) at BJD_TDB=2457170.4990
;                p_3 0.029 -1 # start the fit with Rp/Rstar for planet 3 at 0.020
;                period_3 160 -1 # start the fit with period for planet 3 at 160 days
;              
;                NOTE: if the parameters you put a prior on are not
;                explicitly fit, you must make sure they are
;                appropriately derived in $EXOFAST_PATH/pars2step.pro
;                and $EXOFAST_PATH/step2pars.pro
;
; OPTIONAL INPUTS:
;   NPLANETS    - The number of planets you wish to fit to the
;                 data. Default is 1.
;   RVPATH      - The path to the RV data file(s). The file must have 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) RV (m/s)
;                   3) err (m/s)
;                 NOTE 1: The units must be as specified, or the fit
;                 will be wrong or fail. 
;                 NOTE 2: Other input time stamps will likely not
;                 break the code, but can introduce errors in the
;                 reported times of more than a minute. The output
;                 table will display BJD_TDB but makes no attempt to
;                 convert input times.  See
;                 http://adsabs.harvard.edu/abs/2010PASP..122..935E
;                 for an explanation of times
;                 NOTE 3: If omitted, just the transit data will be
;                 fit.
;   TRANPATH    - The path to the transit data file(s). The file(s) must
;                 have at least 3 columns:
;                   1) Time (BJD_TDB -- See Eastman et al., 2010)
;                   2) Normalized flux
;                   3) err
;                   4) Detrend parameter 1
;                   ....
;                   N+3) Detrend parameter N
;
;                 The names of the files describing the transits
;                 *must* adhere to a certain format:
;                 nYYYYMMDD.filtername.telescope.whateveryouwant. 
;
;                 nYYYYMMDD -- The UTC date of mid transit. This is only
;                 necessary if the data has a single transit. This is
;                 used to label the transits in the output plot.
; 
;                 filtername -- The name of the observed
;                 filter. Only certain values are allowed (use the
;                 closest approximation if yours is not in this list;
;                 see quadld.pro): 
;
;                 Johnson/Cousins: 'U','B','V','R','I','J','H','K'
;                 Sloan: 'Sloanu','Sloang','Sloanr','Sloani','Sloanz'
;                 Kepler: 'Kepler'
;                 CoRoT: 'CoRoT'
;                 Spitzer: 'Spit36','Spit45','Spit58','Spit80'
;                 Stromgren: 'u','b','v','y'
;
;                 This is used to define the limb darkening for the
;                 transit.
;
;                 telescope -- a description of the telescope used for
;                 the observations. Anything is allowed, but all
;                 observations observed with the same telescope should
;                 have the same name. This is used in the output plot
;                 and color codes the TTV plot.
;
;                 whateveryouwant -- any string you want to
;                 include for it to make sense to you. This is not
;                 used by the code.
;
;                 So a transit taken on UTC 2017-01-27 with MINERVA in
;                 the V band would be "n20170127.V.MINERVA.dat"
;
;                 NOTE 2: The units must be as specified, or the fit
;                 will be wrong or fail.
;                 NOTE 3: Other input time stamps will likely not
;                 break the code, but can introduce errors in the
;                 reported times of more than a minute. The output
;                 table will display BJD_TDB but makes no attempt to
;                 convert input times.  See
;                 http://adsabs.harvard.edu/abs/2010PASP..122..935E
;                 for an explanation of times
;                 NOTE 4: If omitted, just the RV data will bit fit
;
;   FLUXFILE   - An ASCII file with each line containing three white
;                space delimited columns: FILTER, MAG,
;                UNCERTAINTY. This file describes the apparent
;                broad-band magnitudes of the host star to fit the SED
;                and derive a distance, luminosity, and
;                extinction. Only certain FILTERS are allowed.  See
;                mag2fluxconv.pro for more info and citations. If not
;                specified, the distance and extinction are not fit
;                and the resultant constraint on the stellar
;                luminosity is not applied. 
;                NOTE: FLUXFILE must be specified for a parallax to
;                constrain the stellar luminosity/radius.
;
;                Galex: 'galNUV', 'galFUV' 
;                NOTE: the SED models are not great at UV, you may find it
;                better not to use Galex
;                WISE: 'WISE1','WISE2','WISE3','WISE4'
;                Sloan: 'uSDSS','gSDSS','rSDSS','iSDSS','zSDSS'
;                Panstarrs: 'zPS'
;                Johnson: 'U','B','V','R',
;                Cousins: 'RC', 'IC','J2M','H2M','K2M'
;                Tycho: 'BT','VT'
;                Kepler INT: 'U_KIS','gKIS','rKIS','iKIS'
;
;                The contents of this file for EPXXXXXXX is below:
; 
;                B      11.949 0.129
;                V      11.572 0.142
;                J2M    10.872 0.022
;                H2M    10.665 0.025
;                K2M    10.601 0.020
;                WISE1  10.537 0.023
;                WISE2  10.583 0.021
;                WISE3  10.618 0.097
;                uSDSS 14.51798 0.05
;                gSDSS 12.14652 0.03
;                rSDSS 11.81735 0.03
;                iSDSS 11.70544 0.03
;                zSDSS 13.20863 0.03
;
;   PREFIX      - Each of the output files will have this string as a
;                 prefix. Default is 'planet.'
;   MAXSTEPS    - The maximum number of steps to take in the MCMC
;                 chain. Note that a 32-bit installation of IDL cannot
;                 allocate more than 260 million double-precision
;                 numbers, and redundant copies of each parameter are
;                 required. Even a 64-bit installation may require
;                 very slow disk swapping. A very large number will
;                 cause memory management problems. Default is
;                 100,000, and larger values are strongly
;                 discouraged. Increase NTHIN if the chains are not
;                 well-mixed. 
;                 NOTE: If the MCMC chains run to MAXSTEPS,
;                 doubling this value will double the runtime. Short
;                 test runs (MAXSTEPS=100) are strongly encouraged
;                 before running week-long fits.
;   NTHIN       - If set, only every NTHINth element will be
;                 kept. Values as high as 1/acceptance rate typically
;                 don't degrade the resultant fit because there is a
;                 high correlation between adjacent steps and has the
;                 advantage of improved memory management and faster
;                 generation of the final plots. Default is 1.
;                 NOTE: Only kept links in the chain count toward
;                 MAXSTEPS, so if the MCMC chains run to MAXSTEPS,
;                 doubling this value will double the runtime.
;   LOGFILE     - A string specifying the name of a file to redirect
;                 the output to. Useful when running many fits or for
;                 a more permanent record of any errors or warnings.
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
;   MAXGR       - The maximum Gelman Rubin statistic that is
;                 considered well-mixed (default=1.01)
;   MINTZ       - The minimum number of independent draws that is
;                 considered well-mixed (default=1000)
;   FITRV       - An NPLANETS boolean array that specifies which
;                 planets should be fit with an RV model. 
;                 If RVPATH is specified, default is bytarr(nplanets) + 1B.
;                 If RVPATH is not specified, default is bytarr(nplanets).
;                 At least one of FITRV and FITTRAN should be true for each planet
;   FITTRAN     - An NPLANETS boolean array that specifies which
;                 planets should be fit with a transit model. 
;                 If TRANPATH is specified, default is bytarr(nplanets) + 1B.
;                 If TRANPATH is not specified, default is bytarr(nplanets). 
;                 At least one of FITRV and FITTRAN should be true for each planet
;   CIRCULAR    - An NPLANETS boolean array that specifies which planets
;                 should be fixed to be circular (1) or left with
;                 eccentricity free (0). 
;                 Default is bytarr(nplanets) + 1B.
;   SECONDARY   - An NPLANETS boolean array that specifies which planets
;                 should have a secondary eclipse modeled. 
;                 Default is bytarr(nplanets)
;                 ***NOT YET IMPLEMENTED***
;   ROSSITER    - An NPLANETS boolean array that specifies which planets
;                 should fit the Rossiter-McLaughlin effect to the RV
;                 data using the Ohta approximation. If a run was
;                 dedicated to RM, it should be separated into its own
;                 file to fit a different zero point and jitter
;                 parameter to it.
;                 Default is bytarr(nplanets)
;                 ***NOT YET IMPLEMENTED***
;   DOPTOM      - An NPLANETS boolean array that specifies which planets
;                 should fit using a Doppler Tomography model.
;                 Default is bytarr(nplanets)
;                 ***NOT YET IMPLEMENTED***
;   CHEN        - An NPLANET boolean array that specifies which
;                 planets should have the Chen & Kipping, 2017
;                 mass-radius relation applied. By default CHEN =
;                 FITRV xor FITTRAN. That is, only apply the
;                 mass-radius prior when RV is not fit (to derive the
;                 planet mass) or when a transit is not fit (to derive
;                 the radius). If the defaults have been overridden
;                 and FITRV and CHEN are both false for a given
;                 planet, the RV semi-amplitude (planet mass) and all
;                 derived parameters will not be quoted. Multi-planet
;                 systems will be constrained not to cross orbits, but
;                 the Hill Sphere will be set to zero.  If FITTRAN and
;                 CHEN are both false for a given planet, the
;                 planetary radius and all derived parameters will not
;                 be quoted.
;
; OPTIONAL KEYWORDS:
;   FITSLOPE  - If set, it will fit a linear trend to the RV data.
;   FITQUAD   - If set, it will fit a quadratic trend to the RV data.
;   NOYY      - If set, do not use the YY evolutionary tracks (see
;               http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2001ApJS..136..417Y)
;               to constrain the mass/radius of the star. If set,
;               \log(rstar) is fit instead of age. This should be set
;               for low-mass stars and an external constraint on the
;               stellar mass and radius should be supplied.
;   TORRES    - If set, use the Torres relations to constrain the mass
;               and radius of the star. This may be useful to
;               investigate potential systematics and should probably
;               be accompanied by the NOYY keyword (but is not
;               enforced). This is not a substitute for YY for low
;               mass stars; the Torres relations are equally
;               inapplicable.
;   NOCLARET  - If set, ignore the Claret & Bloeman limb darkening
;               tables and just fit the limb darkening. This should be
;               specified for low-mass stars where these tables are
;               unreliable.
;   LONGCADENCE - If set, EXPTIME=29.425 and NINTERP=10 are set to handle
;               long cadence data from Kepler. This overwrites
;               NINTERP and EXPTIME inputs.
;   TIDES     - If set, when (1-Rstar/a-rp/a) < e < (1-3*Rstar/a), we
;               set the eccentricity to zero, presuming that the tidal
;               circularization timescale is much much smaller than
;               the age of the system.
;   EPRIOR4   - Parameterize the eccentricity and argument of
;               periastron as e^(1/4)*sin(omega) and
;               e^(1/4)*cos(omega) to more closely match the observed
;               eccentricity distribution
;   TTV       - If set, non-periodic transit times are allowed. The
;               period is constrained by a linear fit to all transit
;               times at each step. Otherwise, a linear ephemeris is
;               assumed.
;               ***NOT YET IMPLEMENTED***
;   TDV       - If set, a new transit depth is fit for each transit
;               ***NOT YET IMPLEMENTED***              
;   BESTONLY  - If set, only the best fit (using AMOEBA) will be
;               performed.
;               ***NOT YET IMPLEMENTED***
;   PLOTONLY  - If set, only a plot of the initial guesses over the
;               data is made. This is useful when refining initial
;               guesses.
;               ***NOT YET IMPLEMENTED***
;   EARTH     - If set, the output units of Mp and Rp are in Earth
;               units, not Jupiter units.
;   DEBUG     - If set, various debugging outputs will be printed and
;               plotted. If the fit is failing, this is your first
;               step to figuring out what's going
;               wrong. Usually a parameter starts too far from its
;               correct value to find or some parameter is not
;               constrained.
; OUTPUTS:
;
;   Each of the output files will be preceeded by PREFIX (defined
;   above).
;
;   mcmc.idl   - An IDL save file that contains the stellar structure
;                with the full chains for each parameter, including
;                derived parameters, the chi2 at each link, and the
;                index of the burn-in period.
;   pdf.ps     - A postscript plot of each posterior distribution
;                function, 8 to a page.
;   covar.ps   - A postscript plot of the covariances between each
;                parameter, 16 to a page. The title of each plot is
;                the correlation coefficient between the two
;                parameters
;   median.tex - The LaTeX source code for a deluxe table of the
;                median values and 68% confidence interval, rounded to
;                two sig figs in the uncertainty.
;   model.ps   - An postscript plot of the best-fit models and residuals.
;                *** Fairly cursory at the moment ***
;   model.rv   - A text file containing the best-fit model RV for each
;                input time
;   model.flux - A text file containing the best-fit model flux for
;                each input time
;
;   NOTE: To extract a single page out of a multi-page PS file, use:
;         psselect -p# input.ps output.ps where # is the page number
;         to extract.
;
; COMMON BLOCKS:
;   CHI2_BLOCK:
;     SS      - A structure that describes the entire stellar system
;
; EXAMPLES:  
;   ;; fit EPXXXXXXX, a 4 planet system with only K2 data
;   exofastv2, nplanets=4, tranpath='ep?????????.Kepler.K2.dat3',fluxfile='ep?????????.flux.txt',$
;           priorfile='ep?????????.priors',debug=debug, prefix='epXXXXXXXXX.eeee.nopilogg.',$
;           fittran=[1,1,1,1],fitrv=[0,0,0,0],circular=[0,0,0,0],/longcadence, /earth
;
; MODIFICATION HISTORY
; 
;  2017/01 -- Complete rewrite of exofast.pro. More general (fits
;             multiple planets, mulitple bands, multiple instrumental
;             offsets). Now easily extensible.
;-
pro exofastv2, priorfile=priorfile, $
               rvpath=rvpath, tranpath=tranpath, fluxfile=fluxfile,$
               prefix=prefix,$
               circular=circular,fitslope=fitslope, secondary=secondary, $
               rossiter=rossiter,$
               nthin=nthin, maxsteps=maxsteps, $
               debug=debug, randomfunc=randomfunc, seed=seed,$
               bestonly=bestonly, plotonly=plotonly,$
               longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
               logfile=logfile, $
               maxgr=maxgr, mintz=mintz, $
               noyy=noyy, tides=tides, nplanets=nplanets, $
               fitrv=fitrv, fittran=fittran,ttv=ttv, earth=earth,$
               i180=i180

;; this is the stellar system structure
COMMON chi2_block, ss

;; name of the chi square function
chi2func = 'exofast_chi2v2'

;; compile all routines now to keep output legible 
;; resolve_all doesn't interpret execute; it's also broken prior to IDL v6.4(?)
if double(!version.release) ge 6.4d0 then $
   resolve_all, resolve_function=[chi2func,'exofast_random'],/cont,/quiet

;; create the master structure
ss = mkss(rvpath=rvpath, tranpath=tranpath, fluxfile=fluxfile, nplanets=nplanets, $
          debug=debug, priorfile=priorfile, fitrv=fitrv, fittran=fittran, $
          circular=circular,fitslope=fitslope, fitquad=fitquad,ttv=ttv, $
          rossiter=rossiter,longcadence=longcadence, earth=earth, i180=i180)

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
print, 'These are the starting values for all the step parameters'
print, 'and the amoeba stepping scale, which is roughly the range'
print, 'of parameter space it will explore around the starting value'
print, 'and is equal to 3x any input priors. When priors are not '
print, 'specified in ' + priorfile + ', a default guess is used.'
print, 'The parameter number is useful for tracking down unconstrained parameters'
print 
print, '**************************************************************'
print, '*** ESPECIALLY DURING BETA TESTING, YOU SHOULD MAKE SURE   ***' 
print, '*** THESE MAKE SENSE AND AGREE WITH YOUR EXPECTATION FROM  ***' 
print, '*** PRIORFILE. IF NOT, YOUR STARTING PRIORS ARE LIKELY NOT ***'
print, '*** BEING TRANSLATED INTO THE STEPPING PARAMETERIZATION    ***'
print, '*** CORRECTLY BY pars2step.pro.                            ***'
print, '**************************************************************'
print
print, 'Par #      Par Name    Par Value       Ameoba Scale'
for i=0, n_elements(name)-1 do print, i, name[i], pars[i], scale[i], format='(i3,x,a15,x,f14.6,x,f14.6)'
;forprint, indgen(n_elements(name)), name, pars, scale,/t,format='(i3,x,a15,x,f14.6,x,f14.6)'
print 

if ss.debug then begin
   print, 'program halted to give you time to inspect the priors. Type ".con" to continue'
   stop
end

print, 'Beginning AMOEBA fit'
best = exofast_amoeba(1d-8,function_name=chi2func,p0=pars,scale=scale,nmax=1d5)
if best[0] eq -1 then begin
   printf, lun, 'ERROR: Could not find best combined fit'
   if lun ne -1 then free_lun, lun
   return
endif
print, 'Finished AMOEBA fit'


;; output the best-fit model fluxes/rvs
bestchi2 = call_function(chi2func,best,modelrv=modelrv,modelflux=modelflux, psname=prefix + 'model.ps')

;; do the MCMC fit
if not keyword_set(bestonly) then begin
   exofast_demc, best, chi2func, pars, chi2=chi2,$
                 nthin=nthin,maxsteps=maxsteps,$
                 burnndx=burnndx, seed=seed, randomfunc=randomfunc, $
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
modelfile = prefix + 'model.mcmc.ps'
bestchi2 = call_function(chi2func,best,psname=modelfile, $
                         modelrv=modelrv, modelflux=modelflux)

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
*(mcmcss.chi2) = chi2

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
