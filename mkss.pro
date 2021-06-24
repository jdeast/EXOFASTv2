;+
; NAME:
;   EXOFAST
;
; PURPOSE:
;
; Creates a stellar system structure that describes an arbitrary
; number of planets, observed bands, rv telescopes (a new zero point
; for each), and observed transits.  This structure is the input to
; basically everything in EXOFASTv2 and is designed to be easily
; extensible. To add a new derived parameter to the output table, add it here,
; then calculate it in DERIVEPARS. To add a new fitted parameter, add
; it here, then use it in EXOFAST_CHI2V2 to influence the model or
; the likelihood/chi^2 directly.
;
; Priors may be applied to any parameter (fitted or derived).
;
; CALLING SEQUENCE:
;   ss = mkss(nplanets=nplanets, circular=circular, $
;             fitslope=fitslope, fitquad=fitquad, ttvs=ttvs, tdvs=tdvs, $
;             rossiter=rossiter, fitdt=fitdt, eprior4=eprior4, fittran=fittran, fitrv=fitrv, $
;             nvalues=nvalues, debug=debug, priorfile=priorfile, $
;             rvpath=rvpath, tranpath=tranpath, longcadence=longcadence, earth=earth)
;
; INPUTS:
; 
;  NPLANETS - The number of planets to fit
;  CIRCULAR - An NPLANETS boolean array describing which
;             planets should be fixed to have circular orbits. If not
;             specified, all planets will be fit with eccentric models.
;  FITTRAN  - An NPLANETS boolean array specifying which planets
;             should be fit with a transit model. By default, all
;             planets are fit with transit photometry. 
;  FITRV    - An NPLANETS boolean array specifying which planets
;             should be fit with a radial velocity model. By default, all
;             planets are fit with radial velocities.
;  CHEN     - An NPLANETS boolean array specifying which planets
;             should have the Chen & Kipping, 2017 mass-radius
;             relation applied. By default CHEN = FITRV xor
;             FITTRAN. That is, only apply the mass-radius prior when
;             RV is not fit (to derive the planet mass) or when a
;             transit is not fit (to derive the radius). If the
;             defaults have been overridden and FITRV and CHEN are
;             both false for a given planet, the RV semi-amplitude
;             (planet mass) and all derived parameters will not be
;             quoted. Multi-planet systems will be constrained not to
;             cross orbits, but the Hill Sphere will be set to zero.
;             If FITTRAN and CHEN are both false for a given planet,
;             the planetary radius and all derived parameters will not
;             be quoted.
;  I180     - An NPLANETS boolean array specifying which planets'
;             inclination should be allowed to be between 0 and 180
;             instead of the default of 0 to 90. Note that a normal
;             transiting planet has a perfect degeneracy between i and
;             180-i which is likely to cause convergence
;             problems. Additional information (e.g., from astrometry
;             or mutual eclipses) must be used for this keyword to be
;             used properly. Even in the case of mutual eclipses
;             (which is currently not supported), at least one planet must be
;             arbitrarily constrained from 0 to 90.
;  ROSSITER - An NPLANETS boolean array specify which planets to fit a
;             rossiter model to the RV data. Fit lambda and
;             V_rot*sin(I_*) from RV data during transit using the
;             Ohta approximations (has known issues)
;  FITDT    - An NPLANETS boolean array specifying which planets to
;             fit a Doppler tomography model to the RV data
;             (*** Currently unsupported***)
;  THERMAL  - A string array specifying which bands to fit thermal
;             emission for. This is what you want to set to fit an
;             isolated secondary eclipse. All observations in this
;             band will be modeled with a baseline of 1 between t2 and
;             t3 of the secondary eclipse and 1 + thermal emission (in
;             PPM) out of eclipse.
;  REFLECT  - A string array specifying which bands to fit reflected
;             light for. Set this along with thermal if you're
;             fitting a full phase curve. It will be modeled as a
;             sinusoid with the orbital period, a minimum at the
;             primary transit, and a fitted amplitude (in PPM).
;  DILUTE   - A string array specifying which bands to fit a dilution
;             term for. Set this if the star is blended with a
;             neighbor and you expect color-dependent depth
;             variations.
;             Note: May be degenerate with F0 (transit normalization) 
;             Note: this only affects the transit model. It is not
;             accounted for in the SED fitting.
;             TODO: automatically model dilution based on multiple
;             SEDs
;  NVALUES  - By default, parameter.value is a scalar. Set this
;             to make it an array of length NVALUES.
;  PRIORFILE - The name of the file that specifies all the priors. The
;              prior file is an ASCII file with each line containing
;              three white space delimited columns: NAME, VALUE,
;              WIDTH. NAME must match a parameter.label. If in an
;              array (e.g., of planets), add "_i", where "i" is the
;              zero-indexed index into the array.  If WIDTH is set to
;              0, the parameter is fixed at VALUE (this is generally
;              not recommended; it's far better to apply a realistic
;              prior). If WIDTH is set to -1, the fit starts at VALUE,
;              but there is no penalty if the model deviates from
;              VALUE. If WIDTH is positive, a gaussian prior is
;              applied. That is, a chi^2 penalty equal to ((parameter
;              - VALUE)/WIDTH)^2 is applied to the likelihood
;              function. Here is the contents of a sample priorfile
;              for the EPXXXXXXXXX system published in Eastman et al,
;              2017:
;
;              teff 6167 78 # Gaussian prior on T_eff of 6167 +/- 78 K
;              feh -0.04 0.1 # Gaussian prior on [Fe/H] of -0.04 +/- 0.1 dex
;              logg 4.22 0.06 # Gaussian prior on logg of 4.22 +/- 0.06
;              vsini 9400 300 # Gaussian prior on vsini of 9400 +/- 300 (ignored since vsini is not fitted)
;              av 0.15 0.30 # Gaussian prior on A_V (extinction) of 0.15 +/- 0.30
;              ma 11.529 0.142 # Gaussian prior on ma (apparent V mag) of 11.529 +/- 0.142
;              distance 247.87 57.223 # Gaussian prior on distance of 247.87 +/- 57.223 pc
;              tc_0 2457166.0543 -1 # start the fit with a TC for planet 0 (EPXXXXXXXXb) at BJD_TDB=2457166.0543
;              p_0 0.020 -1 # start the fit with Rp/Rstar for planet 0 at 0.020
;              period_0 26.847 -1 # start the fit with period for planet 0 at 26.847 days
;              tc_1 2457213.5328 -1 # start the fit with a TC for planet 1 (EPXXXXXXXXc) at BJD_TDB=2457213.5328
;              p_1 0.02 -1 # start the fit with Rp/Rstar for planet 1 at 0.020
;              period_1 39.5419 -1 # start the fit with period for planet 0 at 39.5419 days
;              tc_2 2457191.8894 -1 # start the fit with a TC for planet 2 (EPXXXXXXXXd) at BJD_TDB=2457191.8894
;              p_2 0.020 -1 # start the fit with Rp/Rstar for planet 2 at 0.020
;              period_2 125 -1 # start the fit with period for planet 2 at 125 days
;              tc_3 2457170.4990 -1 # start the fit with a TC for planet 3 (EPXXXXXXXXe) at BJD_TDB=2457170.4990
;              p_3 0.029 -1 # start the fit with Rp/Rstar for planet 3 at 0.020
;              period_3 160 -1 # start the fit with period for planet 3 at 160 days
;  RVPATH   - The path of the RV files. Each RV files is fit with a
;             separate zero point.
;  TRANPATH - The path of the transit data. Each file is fit with a
;             separate normalization. If TTVS is set, each file is fit
;with a diff
;
; KEYWORDS:
;  FITSLOPE - Fit a linear trend to the RV data
;  FITQUAD  - Fit a quadratic trend to the RV data
;  TTVS     - Fit an independent (non-periodic) transit time for each transit
;  TDVS     - Fit an indepdedent depth to each transit
;  EPRIOR4  - Parameterize the eccentricity and argument of
;             periastron as e^(1/4)*sin(omega) and
;             e^(1/4)*cos(omega) to more closely match the observed
;             eccentricity distribution 
;  DEBUG    - Output various debugging information and plots at each
;             step.

; PARAMETER STRUCTURE
; parameter.value -- the parameter's numerical value in
;                    parameter.unit units. For the MCMC, this is an
;                    array for all links.
; parameter.prior -- the parameter's prior
;-
function mkss, nplanets=nplanets, circular=circular,chen=chen, i180=i180,$
               fitslope=fitslope, fitquad=fitquad, $
               ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
               rossiter=rossiter, fitdt=fitdt, eprior4=eprior4, fittran=fittran, fitrv=fitrv, $
               fitdilute=fitdilute, fitthermal=fitthermal, fitreflect=fitreflect, fitphase=fitphase,$
               fitellip=fitellip, fitbeam=fitbeam, derivebeam=derivebeam, $
               nvalues=nvalues, debug=debug, verbose=verbose, priorfile=priorfile, $
               rvpath=rvpath, tranpath=tranpath, dtpath=dtpath, fluxfile=fluxfile, $
               astrompath=astrompath,fitlogmp=fitlogmp,$
               longcadence=longcadence, rejectflatmodel=rejectflatmodel,ninterp=ninterp, exptime=exptime,$
               earth=earth, silent=silent, yy=yy, torres=torres, nomist=nomist, parsec=parsec, $
               noclaret=noclaret,alloworbitcrossing=alloworbitcrossing,$
               logname=logname, best=best, tides=tides, $
               mistsedfile=mistsedfile,fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,oned=oned,$
               fitspline=fitspline, splinespace=splinespace, fitwavelet=fitwavelet, $
               novcve=novcve, nochord=nochord, fitsign=fitsign, randomsign=randomsign, chi2func=chi2func, fittt=fittt,delay=delay, rvepoch=rvepoch

if n_elements(fbolsedfloor) eq 0 then fbolsedfloor = 0.02d0
if n_elements(teffsedfloor) eq 0 then teffsedfloor = 0.015d0
if n_elements(delay) eq 0 then delay = 0L $
else begin
   if delay ne 0 then $
      printandlog, 'You have specified a delay!! This is only designed to test threading and will needlessly slow down the fit. Are you sure?', logname
endelse

if ~keyword_set(nomist) + keyword_set(yy) + keyword_set(parsec) + keyword_set(torres) gt 1 then begin
   printandlog, 'You are **STRONGLY** advised to disable all but one evolutionary model (they are not independent), but type ".con" to proceed', logname
   stop
endif

if n_elements(nplanets) eq 0 then nplanets = 1
if n_elements(circular) ne nplanets and n_elements(circular) gt 1 then begin
   printandlog, "CIRCULAR must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fittran) ne nplanets and n_elements(fittran) gt 1 then begin
   printandlog, "FITTRAN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fitrv) ne nplanets and n_elements(fitrv) gt 1 then begin
   printandlog, "FITRV must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(chen) ne nplanets and n_elements(chen) gt 1 then begin
   printandlog, "CHEN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(i180) ne nplanets and n_elements(i180) gt 1 then begin
   printandlog, "I180 must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(rossiter) ne nplanets and n_elements(rossiter) gt 1 then begin
   printandlog, "ROSSITER must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fitdt) ne nplanets and n_elements(fitdt) gt 1 then begin
   printandlog, "FITDT must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fitlogmp) ne nplanets and n_elements(fitlogmp) gt 1 then begin
   printandlog, "FITLOGMP must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(novcve) ne nplanets and n_elements(novcve) gt 1 then begin
   printandlog, "NOVCVE must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(nochord) ne nplanets and n_elements(nochord) gt 1 then begin
   printandlog, "NOCHORD must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fitsign) ne nplanets and n_elements(fitsign) gt 1 then begin
   printandlog, "FITSIGN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fitbeam) ne nplanets and n_elements(fitbeam) gt 1 then begin
   printandlog, "FITBEAM must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(derivebeam) ne nplanets and n_elements(derivebeam) gt 1 then begin
   printandlog, "DERIVEBEAM must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif
if n_elements(fittt) ne nplanets and n_elements(fittt) gt 1 then begin
   printandlog, "FITTT must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   stop
endif

if not keyword_set(longcadence) then longcadence=0B
if n_elements(fitthermal) eq 0 then fitthermal = ['']
if n_elements(fitreflect) eq 0 then fitreflect = ['']
if n_elements(fitphase) eq 0 then fitphase = ['']
if n_elements(fitellip) eq 0 then fitellip = ['']
if n_elements(fitdilute) eq 0 then fitdilute = ['']
if n_elements(tranpath) eq 0 then tranpath = ''
if n_elements(astrompath) eq 0 then astrompath = ''

;; read in the transit files
if tranpath ne '' or astrompath ne '' then begin

   if tranpath ne '' then begin
      tranfiles=file_search(tranpath,count=ntran)
      if ntran eq 0 then begin
         printandlog, "No transit files files found matching " + strtrim(tranpath,2) + "; please check TRANPATH", logname
         stop
      endif
   endif else begin
      ntran = 0
      tranpath = ''
   endelse

   if astrompath ne '' then begin
      astromfiles=file_search(astrompath,count=nastrom)
      if nastrom eq 0 then begin
         printandlog, "No astrometry files files found matching " + strtrim(astrompath,2) + "; please check ASTROMPATH", logname
         stop
      endif
   endif else begin
      nastrom = 0
      astrompath = ''
   endelse

   ;; find the unique bands
   allowedbands = ['U','B','V','R','I','J','H','K',$
                   'Sloanu','Sloang','Sloanr','Sloani','Sloanz',$
                   'Kepler','TESS','CoRoT','Spit36','Spit45','Spit58','Spit80',$
                   'u','b','v','y']
   prettybands = ['U','B','V','R','I','J','H','K',$
                  "u'","g'","r'","i'","z'",$
                  'Kepler','TESS','CoRoT','$3.6\mu m$','$4.5\mu m$','$5.8\mu m$','$8.0\mu m$',$
                  'u','b','v','y']
   bands = strarr(ntran+nastrom)
   for i=0, ntran+nastrom-1 do begin
      if i lt ntran then bands[i] = (strsplit(file_basename(tranfiles[i]),'.',/extract))(1) $
      else bands[i] = (strsplit(file_basename(astromfiles[i-ntran]),'.',/extract))(1)
      if bands[i] eq 'u' then begin
         bands[i] = 'Sloanu'
      endif else if bands[i] eq 'g' then begin
         bands[i] = 'Sloang'
      endif else if bands[i] eq 'r' then begin
         bands[i] = 'Sloanr'
      endif else if bands[i] eq 'i' then begin
         bands[i] = 'Sloani'
      endif else if bands[i] eq 'z' then begin
         bands[i] = 'Sloanz'
      endif
      if (where(allowedbands eq bands[i]))[0] eq -1 then begin
         printandlog, 'ERROR: band (' + bands[i] + ') not recognized; please select one of the following:'
         printandlog, string(allowedbands)
         stop
      endif
   endfor
   bands = bands[uniq(bands, sort(bands))]
   nband = n_elements(bands)
   if ~keyword_set(silent) then begin
      printandlog, 'The index for each fitted band is', logname
      for i=0L, nband-1 do printandlog, string(i,bands[i],format='(i2,x,a)'),logname
      printandlog, '', logname
   endif
endif else begin
   ntran = 0
   tranpath = ''
   nastrom = 0
   astrompath = ''
   nband = 0
endelse

if nplanets ge 1 and ntran ge 1 then begin
   ;; some error checking on TTVs
   if n_elements(ttvs) eq 0 then ttvs = bytarr(ntran,nplanets) $
   else if n_elements(ttvs) eq 1 then ttvs = bytarr(ntran,nplanets)+ttvs[0] $
   else if n_elements(ttvs) ne ntran*nplanets then begin
      printandlog, 'TTVs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      stop
   endif
   if nplanets gt 1 then begin
      if max(total(ttvs,2)) gt 1 then begin
         printandlog, 'TTVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
   
   ;; some error checking on TIVs
   if n_elements(tivs) eq 0 then tivs = bytarr(ntran,nplanets) $
   else if n_elements(tivs) eq 1 then tivs = bytarr(ntran,nplanets)+tivs[0] $
   else if n_elements(tivs) ne ntran*nplanets then begin
      printandlog, 'TIVs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      stop
   endif
   if nplanets gt 1 then begin
      if max(total(tivs,2)) gt 1 then begin
         printandlog, 'TIVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
   
   ;; some error checking on TDELTAVs
   if n_elements(tdeltavs) eq 0 then tdeltavs = bytarr(ntran,nplanets) $
   else if n_elements(tdeltavs) eq 1 then tdeltavs = bytarr(ntran,nplanets)+tdeltavs[0] $
   else if n_elements(tdeltavs) ne ntran*nplanets then begin
      printandlog, 'Tdeltavs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      stop
   endif
   if nplanets gt 1 then begin
      if max(total(tdeltavs,2)) gt 1 then begin
         printandlog, 'TDELTAVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
endif else begin
   ttvs = 0B
   tivs = 0B
   tdeltavs = 0B
endelse

if n_elements(rejectflatmodel) ne ntran and n_elements(rejectflatmodel) ne 0 then begin
   printandlog, 'REJECTFLATMODEL must be an NTRANSITS element array', logname
   stop
end
if n_elements(rejectflatmodel) eq 0 and ntran gt 0 then rejectflatmodel = bytarr(ntran)

if n_elements(fitspline) ne ntran and n_elements(fitspline) ne 0 then begin
   printandlog, 'FITSPLINE has ' + strtrim(n_elements(fitspline),2) + ' elements; must be an NTRANSITS (' + strtrim(ntran,2) + ') element array', logname
   stop
end
if n_elements(fitspline) eq 0 and ntran gt 0 then fitspline = bytarr(ntran)

if n_elements(splinespace) ne ntran and n_elements(splinespace) ne 0 then begin
   printandlog, 'SPLINESPACE must be an NTRANSITS element array', logname
   stop
end
if n_elements(splinespace) eq 0 and ntran gt 0 then splinespace = dblarr(ntran) + 0.75d0

if n_elements(fitwavelet) ne ntran and n_elements(fitwavelet) ne 0 then begin
   printandlog, 'FITWAVELET must be an NTRANSITS element array', logname
   stop
end
if n_elements(fitwavelet) eq 0 and ntran gt 0 then fitwavelet = bytarr(ntran)


;; if RVPATH not specified or empty, don't use any telescope
if n_elements(rvpath) eq 0 then begin
   ntel = 0
   rvpath = ''
endif else if rvpath eq '' then begin
   ntel = 0
   rvpath = ''
endif else begin
   rvfiles = file_search(rvpath,count=ntel)
endelse

;; was it specifed and not found?
if ntel eq 0 and rvpath ne '' then begin
   printandlog, "RV path (" + rvpath + ") not found! Make sure the file exists or remove the argument to proceed without it." , logname
   stop
endif

;; same for DT path
if n_elements(dtpath) ne 0 then begin
   dtfiles = file_search(dtpath,count=ndt)
   if ndt eq 0 then message, 'DTPATH specified (' + dtpath + ') but no files found'
endif else begin
   dtfiles = ['']
   ndt = 0
endelse

if n_elements(circular) ne nplanets or nplanets eq 0 then circular = bytarr(nplanets>1)
if n_elements(rossiter) ne nplanets or nplanets eq 0 then rossiter = bytarr(nplanets>1)
if n_elements(novcve) ne nplanets or nplanets eq 0 then novcve = bytarr(nplanets>1)
if n_elements(nochord) ne nplanets or nplanets eq 0 then nochord = bytarr(nplanets>1)
if n_elements(fitsign) ne nplanets or nplanets eq 0 then fitsign = bytarr(nplanets>1)+1B
if n_elements(fitbeam) ne nplanets or nplanets eq 0 then fitbeam = bytarr(nplanets>1)
if n_elements(derivebeam) ne nplanets or nplanets eq 0 then derivebeam = bytarr(nplanets>1)

if n_elements(nvalues) ne 0 then value = dblarr(nvalues) $
else value = 0d0
nsteps = n_elements(value)

if n_elements(fittran) eq 0 then begin
   if ntran eq 0 then fittran = bytarr(nplanets>1) $
   else fittran = bytarr(nplanets>1)+1B
endif

if n_elements(fitrv) eq 0 then begin
   if ntel eq 0 then fitrv = bytarr(nplanets>1) $
   else fitrv = bytarr(nplanets>1)+1B
endif

if n_elements(fitdt) eq 0 then begin
   if ndt eq 0 then fitdt = bytarr(nplanets>1) $
   else fitdt = bytarr(nplanets>1)+1B
endif

if n_elements(fittt) eq 0 then begin
   fittt = bytarr(nplanets>1)
endif

if n_elements(fitlogmp) eq 0 then fitlogmp = bytarr(nplanets>1)

if n_elements(chen) ne nplanets or nplanets eq 0 then chen = fittran xor fitrv
if n_elements(i180) ne nplanets or nplanets eq 0 then i180 = bytarr(nplanets>1)

if n_elements(fittran) ne nplanets and nplanets ne 0 then begin
   printandlog, "FITTRAN must have NPLANETS elements", logname
   stop
endif

if n_elements(fitrv) ne nplanets and nplanets ne 0 then begin
   printandlog, "FITRV must have NPLANETS elements", logname
   stop
endif

if nplanets ne 0 then begin
   if (where((~fitrv) and (~fittran)))[0] ne -1 and astrompath eq '' then begin
      printandlog, 'Either a transit or RV must be fit for each planet', logname
      stop
   end
endif

;; was rvpath specified but no planets are fit? (ignore it)
junk = where(fitrv,nrvfit)
if rvpath ne '' and nrvfit eq 0 then begin
   if ~keyword_set(silent) then printandlog, 'WARNING: an RVPATH was specified, but no RV planets are fit. Ignoring the supplied RVs'
   if ~keyword_set(silent) then printandlog, 'Remove RVPATH or set at least one planet in FITRV to 1 to remove this message.'
   ntel = 0
endif

if rvpath eq '' and nrvfit ne 0 then begin
   printandlog, 'ERROR: a RVPATH was not specified, but RVs were requested to be fit (FITRV != 0).'
   printandlog, 'Remove FITRV or specify a RVPATH'
   stop
endif

;; was tranpath specified but no planets are fit? (ignore it)
junk = where(fittran,ntranfit)
if tranpath ne '' and ntranfit eq 0 then begin
   if ~keyword_set(silent) then printandlog, 'WARNING: a TRANPATH was specified, but no transits are fit. Ignoring the supplied transits'
   if ~keyword_set(silent) then printandlog, 'Remove TRANPATH or set at least one planet in FITTRAN to 1 to remove this message.'
   ntran = 0
endif

if tranpath eq '' and ntranfit ne 0 then begin
   printandlog, 'ERROR: a TRANPATH was not specified, but transits were requested to be fit (FITTRAN != 0).'
   printandlog, 'Remove FITTRAN or specify a TRANPATH'
   stop
endif

;; each parameter is a structure, as defined here
parameter = create_struct('value',value,$     ;; its numerical value
                          'prior',0d0,$       ;; its prior value
                          'priorwidth',!values.d_infinity,$ ;; its prior width (infinity => no constraint)
                          'lowerbound',-!values.d_infinity,$ ;; values lower than this have zero likelihood
                          'upperbound',!values.d_infinity,$ ;; values higher than this have zero likelihood
                          'label','',$        ;; what do I call it?
                          'cgs',1d0,$         ;; multiply value by this to convert to cgs units
                          'link',ptr_new(),$  ;; a pointer to a linked parameter (not implemented)
                          'scale',0d0,$       ;; scale for amoeba
                          'latex','',$        ;; latex label for the latex table
                          'userchanged',0B,$  ;; the user supplied a prior that impacts this parameter
                          'description','',$  ;; a verbose description
                          'unit','',$         ;; units ("Radians" triggers special angular behavior)
                          'medvalue','',$     ;; median value
                          'upper','',$        ;; upper error bar (68% confidence)
                          'lower','',$        ;; lower error bar (68% confidence)
                          'best',!values.d_nan,$ ;; best-fit parameter
                          'fit', 0B,$         ;; If true, step in this parameter during the MCMC
                          'derive',1B)        ;; If true, quote this parameter in the final table
                          
;; define each parameter (derived, fitted, and ignored)
logmstar = parameter
logmstar.value = 0d0
logmstar.unit = ''
logmstar.description = 'log Mass'
logmstar.latex = '\log{M_*/M_{\sun}}'
logmstar.label = 'logmstar'
logmstar.cgs = !values.d_nan
logmstar.fit = 1
logmstar.derive=0
logmstar.scale = 1

mstar = parameter
mstar.value = 10^logmstar.value
mstar.unit = '\msun'
mstar.description = 'Mass'
mstar.latex = 'M_*'
mstar.label = 'mstar'
mstar.cgs = 1.9891d33

logg = parameter
logg.value = 4.44d0
logg.unit = 'cgs'
logg.description = 'Surface gravity'
logg.latex = '\log{g}'
logg.label = 'logg'
logg.scale = 0.3d0

teff = parameter
teff.value = 5778d0
teff.unit = 'K'
teff.description = 'Effective Temperature'
teff.latex = 'T_{\rm eff}'
teff.label = 'teff'
teff.fit=1
teff.scale = 500d0

teffsed = parameter
teffsed.value = 5778d0
teffsed.unit = 'K'
teffsed.description = 'Effective Temperature'
teffsed.latex = 'T_{\rm eff,SED}'
teffsed.label = 'teffsed'
teffsed.fit=0
teffsed.derive=0
teffsed.scale = 500d0

feh = parameter
feh.value = 0d0
feh.description = 'Metallicity'
feh.latex = '[{\rm Fe/H}]'
feh.label = 'feh'
feh.unit = 'dex'
feh.fit=1
feh.scale = 0.5d0

initfeh = parameter
initfeh.value = 0d0
initfeh.description = 'Initial Metallicity'
initfeh.latex = '[{\rm Fe/H}]_{0}'
initfeh.label = 'initfeh'
initfeh.scale = 0.5d0
if ~keyword_set(nomist) or keyword_set(parsec) then begin
   initfeh.fit=1
   initfeh.derive=1
endif else begin
   initfeh.fit=0
   initfeh.derive=0
endelse

Av = parameter
Av.description = 'V-band extinction'
Av.latex = 'A_V'
Av.label = 'Av'
Av.unit = 'mag'
Av.fit = 0
Av.derive = 0
Av.scale = 0.3d0

alpha = parameter
alpha.description = 'Alpha abundance'
alpha.latex = '\alpha'
alpha.label = 'alpha'
alpha.fit = 0
alpha.derive = 0
alpha.scale = 0.3d0

;Ma = parameter
;Ma.unit = ''
;Ma.description = 'Apparent V-band Magnitude'
;Ma.latex = 'M_a'
;Ma.label = 'Ma'
;Ma.fit=1

;Mv = parameter
;Mv.unit = ''
;Mv.description = 'Absolute V-band Magnitude'
;Mv.latex = 'M_v'
;Mv.label = 'Mv'

distance = parameter
distance.unit = 'pc'
distance.description = 'Distance'
distance.latex = 'd'
distance.label = 'distance'
distance.cgs = 3.08567758d18 ;; cm/pc
distance.scale = 100
distance.value = 10
distance.fit = 0B
distance.derive = 0B
if nastrom gt 0 then begin
   distance.fit = 1d0
   distance.derive = 1B
endif

fbol = parameter
fbol.unit = 'cgs'
fbol.description = 'Bolometric Flux'
fbol.latex = 'F_{Bol}'
fbol.label = 'fbol'
fbol.cgs = 1d0
fbol.derive = 0d0

ra = parameter
ra.value = 0d0
ra.scale = 0.1d0
ra.unit = 'ICRS deg'
ra.description = 'Right Ascension'
ra.latex = '\alpha_{ICRS}'
ra.label = 'ra'
ra.cgs = 3600d0*180d0/!dpi ;; rad/deg
ra.derive = 0B

dec = parameter
dec.value = 0d0
dec.scale = 0.1d0
dec.unit = 'ICRS deg'
dec.description = 'Declination'
dec.latex = '\delta_{ICRS}'
dec.label = 'dec'
dec.cgs = 3600d0*180d0/!dpi ;; rad/deg
dec.derive = 0B

pmra = parameter
pmra.value = 0d0
pmra.scale = 1000d0
pmra.unit = 'mas/yr'
pmra.description = 'RA Proper Motion'
pmra.latex = '\mu_{\alpha*}'
pmra.label = 'pmra'
pmra.cgs = 3600d3*180d0/!dpi*86400d0*365.25d0 ;; rad/deg
pmra.derive = 0B

pmdec = parameter
pmdec.value = 0d0
pmdec.scale = 1000d0
pmdec.unit = 'mas/yr'
pmdec.description = 'Dec Proper Motion'
pmdec.latex = '\mu_{\delta}'
pmdec.label = 'pmdec'
pmdec.cgs = 3600d3*180d0/!dpi*86400d0*365.25d0 ;; rad/deg
pmdec.derive = 0B

rvabs = parameter
rvabs.value = 0d0
rvabs.scale = 10000d0
rvabs.unit = 'm/s'
rvabs.description = 'Absolute RV'
rvabs.latex = '\gamma_{\rm abs}'
rvabs.label = 'rvabs'
rvabs.cgs = 100d0 ;; m/s
rvabs.derive = 0B

astromscale = parameter
astromscale.value = 1d0
astromscale.scale = 10d0
astromscale.unit = ''
astromscale.description = 'Astrometric Error Scaling'
astromscale.latex = '\sigma_{Astrom}'
astromscale.label = 'astromscale'
astromscale.cgs = 0
if nastrom gt 0 then begin
   astromscale.fit = 1B
   astromscale.derive=1B
endif else astromscale.derive = 0B

raoffset = parameter
raoffset.value = 0d0
raoffset.scale = 1d-2
raoffset.unit = 'degrees'
raoffset.description = 'RA offset'
raoffset.latex = '\Delta\alpha'
raoffset.label = 'raoffset'
raoffset.cgs = 0 ;;
raoffset.derive = 0B

decoffset = parameter
decoffset.value = 0d0
decoffset.scale = 1d-2
decoffset.unit = 'degrees'
decoffset.description = 'Dec offset'
decoffset.latex = '\Delta\delta'
decoffset.label = 'decoffset'
decoffset.cgs = 1d0 ;;
decoffset.derive = 0B

phottobary = parameter
phottobary.value = 1d0
phottobary.scale = 10d0
phottobary.unit = ''
phottobary.description = 'Photocenter/Barycenter'
phottobary.latex = 'p/\rho'
phottobary.label = 'phottobary'
phottobary.cgs = 1d0
phottobary.derive = 0B

parallax = parameter
parallax.value = 0d0
parallax.scale = 100d0
parallax.unit = 'mas'
parallax.description = 'Parallax'
parallax.latex = '\varpi'
parallax.label = 'parallax'
parallax.cgs = 3600d3*180d0/!dpi ;; rad/mas
parallax.derive=0B
parallax.fit = 0B
if nastrom gt 0 then begin
   parallax.derive=1B
   fbol.derive = 1B
endif

gamma = parameter
gamma.unit = 'm/s'
gamma.description = 'Relative RV Offset'
gamma.latex = '\gamma_{\rm rel}'
gamma.label = 'gamma'
gamma.cgs = 100d0
if ntel gt 0 and nplanets gt 0 then gamma.fit=1 $
else gamma.derive=0
gamma.scale=5000d0

slope = parameter
slope.unit = 'm/s/day'
slope.description = 'RV slope'
slope.latex = '\dot{\gamma}'
slope.label = 'slope'
slope.cgs = 100d0/86400d0
if keyword_set(fitslope) or keyword_set(fitquad) then slope.fit = 1 $
else slope.derive=0
slope.scale=1d0

quad = parameter
quad.unit = 'm/s/day^2'
quad.description = 'RV quadratic term'
quad.latex = '\ddot{\gamma}'
quad.label = 'quad'
quad.cgs = 100d0/86400d0^2
if keyword_set(fitquad) then quad.fit=1 $
else quad.derive=0
quad.scale = 1d0

rstar = parameter
rstar.value = 1d0
rstar.unit = '\rsun'
rstar.description = 'Radius'
rstar.latex = 'R_*'
rstar.label = 'rstar'
rstar.cgs = 6.955d10
rstar.fit = 1
rstar.scale = 0.5d0

rstarsed = parameter
rstarsed.value = 1d0
rstarsed.unit = '\rsun'
rstarsed.description = 'Radius'
rstarsed.latex = 'R_{*,SED}'
rstarsed.label = 'rstarsed'
rstarsed.cgs = 6.955d10
rstarsed.fit = 0
rstarsed.derive = 0
rstarsed.scale = 0.5d0

age = parameter
age.value = 4.603d0
age.unit = 'Gyr'
age.description = 'Age'
age.latex = 'Age'
age.label = 'age'
age.cgs = 3600d0*24d0*365.242d0*1d9
age.derive = 0
if keyword_set(yy) then begin
   age.fit = 1
   age.derive = 1
endif
if ~keyword_set(nomist) or keyword_set(parsec) then age.derive=1

age.scale = 3d0

eep = parameter
eep.value = 398.668d0 ;; solar EEP (PARSEC)
eep.value = 354.1661d0 ;; solar EEP (MIST)
eep.unit = ''
eep.description = 'Equal Evolutionary Phase'
eep.latex = 'EEP'
eep.label = 'eep'
eep.scale = 50d0
if ~keyword_set(nomist) or keyword_set(parsec) then begin
   eep.fit = 1
   eep.derive = 1
endif else begin
   eep.derive=0
   eep.fit=0
endelse

lstar = parameter
lstar.unit = '\lsun'
lstar.description = 'Luminosity'
lstar.latex = 'L_*'
lstar.label = 'lstar'
lstar.cgs = 3600d0*24d0*365.242d0*1d9

rhostar = parameter
rhostar.unit = 'cgs'
rhostar.description = 'Density'
rhostar.latex = '\rho_*'
rhostar.label = 'rhostar'

tc = parameter
tc.unit = '\bjdtdb'
tc.description = 'Time of conjunction'
tc.latex = 'T_C'
tc.label = 'tc'
tc.cgs = 86400d0
if nplanets eq 0 then tc.derive=0 $
else tc.fit = 1
tc.scale = 0.1

tt = parameter
tt.unit = '\bjdtdb'
tt.description = 'Time of minimum projected separation'
tt.latex = 'T_T'
tt.label = 'tt'
tt.cgs = 86400d0
if ntran eq 0 then tt.derive=0
tt.scale = 0.1
tt.derive = 0B

t0 = parameter
t0.unit = '\bjdtdb'
t0.description = 'Optimal conjunction Time'
t0.latex = 'T_0'
t0.label = 't0'
t0.cgs = 86400d0
if nplanets eq 0 then t0.derive=0

logp = parameter
logp.unit = ''
logp.description = 'Log of Period'
logp.latex = '\log{P}'
logp.label = 'logp'
logp.cgs = !values.d_nan
if nplanets ne 0 then logp.fit = 1
logp.derive=0
logp.scale = 0.01d0

qesinw = parameter
qesinw.unit = ''
qesinw.description = ''
qesinw.latex = 'e^{1/4}\sin{\omega_*}'
qesinw.label = 'qesinw'
qesinw.cgs = !values.d_nan
qesinw.scale = 0.1d0
qesinw.derive = 0

qecosw = parameter
qecosw.unit = ''
qecosw.description = ''
qecosw.latex = 'e^{1/4}\cos{\omega_*}'
qecosw.label = 'qecosw'
qecosw.cgs = !values.d_nan
qecosw.scale = 0.1d0
qecosw.derive = 0

sesinw = parameter
sesinw.unit = ''
sesinw.description = ''
sesinw.latex = '\sqrt{e}\sin{\omega_*}'
sesinw.label = 'sesinw'
sesinw.cgs = !values.d_nan
sesinw.scale = 0.1d0
sesinw.derive = 0

secosw = parameter
secosw.unit = ''
secosw.description = ''
secosw.latex = '\sqrt{e}\cos{\omega_*}'
secosw.label = 'secosw'
secosw.cgs = !values.d_nan
secosw.scale = 0.1d0
secosw.derive = 0

vcve = parameter
vcve.unit = ''
vcve.description = ''
vcve.latex = 'V_c/V_e'
vcve.label = 'vcve'
vcve.cgs = !values.d_nan
vcve.scale = 1d0
vcve.value = 0.99d0
if nplanets eq 0 then vcve.derive=0B

chord = parameter
chord.unit = ''
chord.description = ''
chord.latex = '((1-R_P/R_*)^2-b^2)^{1/2}'
chord.label = 'chord'
chord.cgs = 1d0
chord.scale = 1d0
chord.value = 1d0
chord.derive = 0

sign = parameter
sign.unit = ''
sign.description = ''
sign.latex = 'sign'
sign.label = 'sign'
sign.cgs = !values.d_nan
sign.scale = 1
sign.value = 1.5
sign.derive = 0

sign2 = parameter
sign2.unit = ''
sign2.description = ''
sign2.latex = 'sign2'
sign2.label = 'sign2'
sign2.cgs = !values.d_nan
sign2.scale = 1
sign2.value = 1.5
sign2.derive = 0

if keyword_set(eprior4) then begin
   ;; step in qesinw (prior favoring lower e)
   qesinw.fit = 1
   qecosw.fit = 1
endif else begin
   ;; step in sqrt(e)sinw (uniform prior)
   if nplanets ne 0 then begin
      sesinw.fit = 1
      secosw.fit = 1
   endif
endelse

esinw = parameter
esinw.unit = ''
esinw.description = ''
esinw.latex = 'e\sin{\omega_*}'
esinw.label = 'esinw'
esinw.cgs = !values.d_nan
if nplanets eq 0 then esinw.derive=0

ecosw = parameter
ecosw.unit = ''
ecosw.description = ''
ecosw.latex = 'e\cos{\omega_*}'
ecosw.label = 'ecosw'
ecosw.cgs = !values.d_nan
if nplanets eq 0 then ecosw.derive=0

e = parameter
e.unit = ''
e.description = 'Eccentricity'
e.latex = 'e'
e.label = 'e'
e.cgs = !values.d_nan
if nplanets eq 0 then e.derive=0

omega = parameter
omega.unit = 'Radians'
omega.description = 'Argument of Periastron'
omega.latex = '\omega_*'
omega.label = 'omega'
omega.cgs = 1d0
omega.derive = 0

lsinw = parameter
lsinw.unit = ''
lsinw.description = 'L*Sine of Argument of Periastron'
lsinw.latex = 'L\sin\{\omega_*}'
lsinw.label = 'lsinw'
lsinw.cgs = 1d0
lsinw.value = 0d0
lsinw.scale = 1d0
lsinw.derive = 0

lsinw2 = parameter
lsinw2.unit = ''
lsinw2.description = 'L*Sine of Argument of Periastron'
lsinw2.latex = 'sign(vcve-1)*L\sin\{\omega_*}'
lsinw2.label = 'lsinw2'
lsinw2.cgs = 1d0
lsinw2.value = 0d0
lsinw2.scale = 1d0
lsinw2.derive = 0

lcosw = parameter
lcosw.unit = ''
lcosw.description = 'L*Cosine of Argument of Periastron'
lcosw.latex = 'L\cos{\omega_*}'
lcosw.label = 'lcosw'
lcosw.cgs = 1d0
lcosw.value = 0d0
lcosw.scale = 1d0
lcosw.derive = 0

omegadeg = parameter
omegadeg.unit = 'Degrees'
omegadeg.description = 'Argument of Periastron'
omegadeg.latex = '\omega_*'
omegadeg.label = 'omegadeg'
omegadeg.cgs = !dpi/180d0
if nplanets eq 0 then omegadeg.derive=0

bigomega = parameter
bigomega.value = 0d0
bigomega.scale = !dpi
bigomega.unit = 'Radians'
bigomega.description = 'Longitude of ascending node'
bigomega.latex = '\Omega'
bigomega.label = 'bigomega'
bigomega.cgs = 1d0
bigomega.derive = 0B

lsinbigomega = parameter
lsinbigomega.value = 0d0
lsinbigomega.scale = 0.5
lsinbigomega.unit = ''
lsinbigomega.description = 'l*Sine of Longitude of ascending node'
lsinbigomega.latex = 'l\sin{\Omega}'
lsinbigomega.label = 'lsinbigomega'
lsinbigomega.cgs = 1d0
lsinbigomega.derive = 0B

lcosbigomega = parameter
lcosbigomega.value = 0d0
lcosbigomega.scale = 0.5
lcosbigomega.unit = ''
lcosbigomega.description = 'l*Cosine of Longitude of ascending node'
lcosbigomega.latex = 'l\cos{\Omega}'
lcosbigomega.label = 'lcosbigomega'
lcosbigomega.cgs = 1d0
lcosbigomega.derive = 0B

bigomegadeg = parameter
bigomegadeg.value = 0d0
bigomegadeg.scale = 180d0
bigomegadeg.unit = 'Degrees'
bigomegadeg.description = 'Longitude of ascending node'
bigomegadeg.latex = '\Omega'
bigomegadeg.label = 'bigomegadeg'
bigomegadeg.cgs = !dpi/180d0
bigomegadeg.derive=0B
bigomegadeg.fit=0B
if nastrom gt 0 then begin
   lsinbigomega.fit=1B
   lcosbigomega.fit=1B
   bigomegadeg.derive=1B
endif else begin
   lsinbigomega.fit=0B
   lcosbigomega.fit=0B
endelse
   
p = parameter
p.value = 0.1d0
p.unit = ''
p.description = 'Radius of planet in stellar radii'
p.latex = 'R_P/R_*'
p.label = 'p'
p.scale = 1d-1
if nplanets eq 0 then p.derive = 0 $
else p.fit = 1

ar = parameter
ar.unit = ''
ar.description = 'Semi-major axis in stellar radii'
ar.latex = 'a/R_*'
ar.label = 'ar'
if nplanets eq 0 then ar.derive=0

cosi = parameter
cosi.value = 0d0
cosi.unit = ''
cosi.description = 'Cos of inclination'
cosi.latex = '\cos{i}'
cosi.label = 'cosi'
if nplanets ne 0 then cosi.fit = 1
cosi.scale = 0.1
cosi.derive = 0

i = parameter
i.unit = 'Radians'
i.description = 'Inclination'
i.latex = 'i'
i.label = 'i'
i.derive = 0

ideg = parameter
ideg.unit = 'Degrees'
ideg.description = 'Inclination'
ideg.latex = 'i'
ideg.label = 'ideg'
ideg.cgs = !dpi/180d0
if nplanets eq 0 then ideg.derive = 0

tp = parameter
tp.unit = '\bjdtdb'
tp.description = 'Time of Periastron'
tp.latex = 'T_P'
tp.label = 'tp'
tp.cgs = 86400d0
if nplanets eq 0 then tp.derive = 0

td = parameter
td.unit = '\bjdtdb'
td.description = 'Time of Descending Node'
td.latex = 'T_D'
td.label = 'td'
td.cgs = 86400d0
if nplanets eq 0 then td.derive = 0

ta = parameter
ta.unit = '\bjdtdb'
ta.description = 'Time of Ascending Node'
ta.latex = 'T_A'
ta.label = 'ta'
ta.cgs = 86400d0
if nplanets eq 0 then ta.derive = 0

ts = parameter
ts.unit = '\bjdtdb'
ts.description = 'Time of eclipse'
ts.latex = 'T_S'
ts.label = 'ts'
ts.cgs = 86400d0
if nplanets eq 0 then ts.derive = 0

phase = parameter
phase.unit = ''
phase.description = 'Phase of inferior conjunction (primary transit)'
phase.latex = '\phi'
phase.label = 'phase'
phase.derive = 0

tfwhm = parameter
tfwhm.unit = 'days'
tfwhm.description = 'FWHM transit duration'
tfwhm.latex = 'T_{FWHM}'
tfwhm.label = 'tfwhm'
tfwhm.cgs = 86400d0
if nplanets eq 0 then tfwhm.derive = 0

t14 = parameter
t14.unit = 'days'
t14.description = 'Total transit duration'
t14.latex = 'T_{14}'
t14.label = 't14'
t14.cgs = 86400d0
if nplanets eq 0 then t14.derive = 0

tau = parameter
tau.unit = 'days'
tau.description = 'Ingress/egress transit duration'
tau.latex = '\tau'
tau.label = 'tau'
tau.cgs = 86400d0
if nplanets eq 0 then tau.derive = 0

tfwhms = parameter
tfwhms.unit = 'days'
tfwhms.description = 'FWHM eclipse duration'
tfwhms.latex = 'T_{S,FWHM}'
tfwhms.label = 'tfwhms'
tfwhms.cgs = 86400d0
if nplanets eq 0 then tfwhms.derive = 0

t14s = parameter
t14s.unit = 'days'
t14s.description = 'Total eclipse duration'
t14s.latex = 'T_{S,14}'
t14s.label = 't14s'
t14s.cgs = 86400d0
if nplanets eq 0 then t14s.derive = 0

taus = parameter
taus.unit = 'days'
taus.description = 'Ingress/egress eclipse duration'
taus.latex = '\tau_S'
taus.label = 'taus'
taus.cgs = 86400d0
if nplanets eq 0 then taus.derive = 0

eclipsedepth = parameter
eclipsedepth.unit = 'ppm'
eclipsedepth.description = 'Measured eclipse depth'
eclipsedepth.latex = '\delta_{S}'
eclipsedepth.label = 'eclipsedepth'
eclipsedepth.cgs = 1d6
eclipsedepth.derive = 0

;eclipsedepth36 = parameter
;eclipsedepth36.unit = 'ppm'
;eclipsedepth36.description = 'Blackbody eclipse depth at 3.6$\mu$m'
;eclipsedepth36.latex = '\delta_{S,3.6\mu m}'
;eclipsedepth36.label = 'eclipsedepth36'
;eclipsedepth36.cgs = 1d6
;if nplanets eq 0 then eclipsedepth36.derive = 0

;eclipsedepth45 = parameter
;eclipsedepth45.unit = 'ppm'
;eclipsedepth45.description = 'Blackbody eclipse depth at 4.5$\mu$m'
;eclipsedepth45.latex = '\delta_{S,4.5\mu m}'
;eclipsedepth45.label = 'eclipsedepth45'
;eclipsedepth45.cgs = 1d6
;if nplanets eq 0 then eclipsedepth45.derive = 0

eclipsedepth25 = parameter
eclipsedepth25.unit = 'ppm'
eclipsedepth25.description = 'Blackbody eclipse depth at 2.5$\mu$m'
eclipsedepth25.latex = '\delta_{S,2.5\mu m}'
eclipsedepth25.label = 'eclipsedepth25'
eclipsedepth25.cgs = 1d6
if nplanets eq 0 then eclipsedepth25.derive = 0

eclipsedepth50 = parameter
eclipsedepth50.unit = 'ppm'
eclipsedepth50.description = 'Blackbody eclipse depth at 5.0$\mu$m'
eclipsedepth50.latex = '\delta_{S,5.0\mu m}'
eclipsedepth50.label = 'eclipsedepth50'
eclipsedepth50.cgs = 1d6
if nplanets eq 0 then eclipsedepth50.derive = 0

eclipsedepth75 = parameter
eclipsedepth75.unit = 'ppm'
eclipsedepth75.description = 'Blackbody eclipse depth at 7.5$\mu$m'
eclipsedepth75.latex = '\delta_{S,7.5\mu m}'
eclipsedepth75.label = 'eclipsedepth75'
eclipsedepth75.cgs = 1d6
if nplanets eq 0 then eclipsedepth75.derive = 0

delta = parameter
delta.unit = 'fraction'
delta.description = 'Transit depth'
delta.latex = '\delta'
delta.label = 'delta'
if nplanets eq 0 then delta.derive = 0

depth = parameter
depth.description = 'Flux decrement at mid transit'
depth.latex = 'Depth'
depth.label = 'depth'
if nplanets eq 0 then depth.derive = 0

dr = parameter
dr.description = 'Separation at mid transit'
dr.latex = 'd/R_*'
dr.label = 'dr'
if nplanets eq 0 then dr.derive = 0

lambda = parameter
lambda.unit = 'Radians'
lambda.description = 'Projected Spin-orbit alignment'
lambda.latex = '\lambda'
lambda.label = 'lambda'
lambda.scale = !dpi
lambda.derive = 0

lsinlambda = parameter
lsinlambda.value = 0d0
lsinlambda.unit = ''
lsinlambda.description = 'lSine of Projected Spin-orbit alignment'
lsinlambda.latex = 'l\sin{\lambda}'
lsinlambda.label = 'lsinlambda'
lsinlambda.scale = 0.5
lsinlambda.derive = 0

lcoslambda = parameter
lcoslambda.value = 0d0
lcoslambda.unit = ''
lcoslambda.description = 'l*Cosine of Projected Spin-orbit alignment'
lcoslambda.latex = 'l\cos{\lambda}'
lcoslambda.label = 'lcoslambda'
lcoslambda.scale = 0.5
lcoslambda.derive = 0

lambdadeg = parameter
lambdadeg.unit = 'Degrees'
lambdadeg.description = 'Projected Spin-orbit alignment'
lambdadeg.latex = '\lambda'
lambdadeg.label = 'lambdadeg'
lambdadeg.cgs = !dpi/180d0
lambdadeg.derive = 0

vsini = parameter
vsini.unit = 'm/s'
vsini.description = 'Projected rotational velocity'
vsini.latex = 'vsinI_*'
vsini.label = 'vsini'
vsini.cgs = 100d0
vsini.derive = 0
vsini.scale = 5d3

vline = parameter
vline.unit = 'm/s'
vline.description = 'Unbroadened line width'
vline.latex = 'V_{\rm line}'
vline.label = 'vline'
vline.cgs = 1000d0
vline.derive = 0
vline.value = 4d3
vline.scale = 2d3

dtscale = parameter
dtscale.description = 'Doppler Tomography Error scaling'
dtscale.latex = '\sigma_{DT}'
dtscale.label = 'dtscale'
dtscale.derive = 0
dtscale.scale = 1d2
dtscale.value = 1d0

logmp = parameter
logmp.description = 'Log of mp in solar units'
logmp.latex = '\log{M_P}'
logmp.label = 'logmp'
logmp.derive=0

k = parameter
k.value = 10d0
k.unit = 'm/s'
k.description = 'RV semi-amplitude'
k.latex = 'K'
k.label = 'k'
k.cgs = 100d0
if nplanets eq 0 then k.derive=0

logk = parameter
logk.value = 10d0
logk.unit = 'm/s'
logk.description = 'RV semi-amplitude'
logk.latex = 'LOGK'
logk.label = 'logk'
logk.cgs = 100d0
logk.derive=0

period = parameter
period.unit = 'days'
period.description = 'Period'
period.latex = 'P'
period.label = 'Period'
if nplanets eq 0 then period.derive=0

a = parameter
a.unit = 'AU'
a.description = 'Semi-major axis'
a.latex = 'a'
a.label = 'a'
if nplanets eq 0 then a.derive=0

arsun = parameter
arsun.unit = '\rsun'
arsun.description = 'Semi-major axis in solar radii'
arsun.latex = 'arsun'
arsun.label = 'arsun'
arsun.cgs = 6.955d10
arsun.derive = 0

f0 = parameter
f0.value = 1d0
f0.description = 'Baseline flux'
f0.latex = 'F_0'
f0.label = 'f0'
f0.scale = 1d-2
if nplanets eq 0 then f0.derive = 0 $
else f0.fit = 1

detrend = parameter
detrend.value = 0d0
detrend.description = 'Detrending variable'
detrend.latex = ''
detrend.label = 'detrend'
detrend.scale = 0.1
detrend.fit = 1
detrend.derive = 1

beam = parameter
beam.value = 0d0
beam.description = 'Doppler Beaming amplitude'
beam.latex = 'A_B'
beam.unit = 'ppm'
beam.label = 'beam'
beam.scale = 100d0
beam.fit = 0
beam.derive = 0

bs = parameter
bs.description = 'Eclipse impact parameter'
bs.latex = 'b_S'
bs.label = 'bs'
if nplanets eq 0 then bs.derive = 0

ps = parameter
ps.description = 'A priori non-grazing eclipse prob'
ps.latex = 'P_S'
ps.label = 'ps'
if nplanets eq 0 then ps.derive = 0

psg = parameter
psg.description = 'A priori eclipse prob'
psg.latex = 'P_{S,G}'
psg.label = 'psg'
if nplanets eq 0 then psg.derive = 0

pt = parameter
pt.description = 'A priori non-grazing transit prob'
pt.latex = 'P_T'
pt.label = 'pt'
if nplanets eq 0 then pt.derive = 0

ptg = parameter
ptg.description = 'A priori transit prob'
ptg.latex = 'P_{T,G}'
ptg.label = 'ptg'
if nplanets eq 0 then ptg.derive = 0

mp = parameter
mp.unit = '\mj'
mp.description = 'Mass'
mp.latex = 'M_P'
mp.label = 'mp'
mp.cgs = 1.89813d30
if nplanets eq 0 then mp.derive = 0

mpearth = parameter
mpearth.unit = '\me'
mpearth.description = 'Mass'
mpearth.latex = 'M_P'
mpearth.label = 'mpearth'
mpearth.cgs = 5.97219d27
if nplanets eq 0 then mpearth.derive = 0

mpsun = parameter
mpsun.unit = '\msun'
mpsun.description = 'Mass'
mpsun.latex = 'M_P'
mpsun.label = 'mpsun'
mpsun.cgs = 1.9891d33
mpsun.derive = 0
mpsun.value = 0.001d0
if nplanets ne 0 then mpsun.fit = 1

msini = parameter
msini.unit = '\mj'
msini.description = 'Minimum mass'
msini.latex = 'M_P\sin i'
msini.label = 'msini'
msini.cgs = 1.89813d30
if nplanets eq 0 then msini.derive = 0

msiniearth = parameter
msiniearth.unit = '\me'
msiniearth.description = 'Minimum mass'
msiniearth.latex = 'M_P\sin i'
msiniearth.label = 'msiniearth'
msiniearth.cgs = 5.97219d27
if nplanets eq 0 then msiniearth.derive = 0

b = parameter
b.description = 'Transit Impact parameter'
b.latex = 'b'
b.label = 'b'
if nplanets eq 0 then b.derive = 0

q = parameter
q.description = 'Mass ratio'
q.latex = 'M_P/M_*'
q.label = 'q'
if nplanets eq 0 then q.derive = 0

rp = parameter
rp.unit = '\rj'
rp.description = 'Radius'
rp.latex = 'R_P'
rp.label = 'rp'
rp.cgs = 7.1492d9
if nplanets eq 0 then rp.derive = 0

rpsun = parameter
rpsun.unit = '\rsun'
rpsun.description = 'Radius'
rpsun.latex = 'R_P'
rpsun.label = 'rpsun'
rpsun.cgs = 6.955d10
rpsun.derive = 0

rpearth = parameter
rpearth.unit = '\re'
rpearth.description = 'Radius'
rpearth.latex = 'R_P'
rpearth.label = 'rpearth'
rpearth.cgs = 6.3781d8
if nplanets eq 0 then rpearth.derive = 0

if keyword_set(earth) then begin
   ;; don't display jupiter unit parameters
   rp.derive = 0
   mp.derive = 0
   msini.derive = 0
endif else begin
   ;; don't display earth unit parameters
   rpearth.derive = 0
   mpearth.derive = 0
   msiniearth.derive = 0
endelse

rhop = parameter
rhop.unit = 'cgs'
rhop.description = 'Density'
rhop.latex = '\rho_P'
rhop.label = 'rhop'
if nplanets eq 0 then rhop.derive = 0

loggp = parameter
loggp.description = 'Surface gravity'
loggp.latex = 'logg_P'
loggp.label = 'loggp'
if nplanets eq 0 then loggp.derive = 0

teq = parameter
teq.unit = 'K'
teq.description = 'Equilibrium temperature'
teq.latex = 'T_{eq}'
teq.label = 'teq'
if nplanets eq 0 then teq.derive = 0

tcirc = parameter
tcirc.unit = 'Gyr'
tcirc.description = 'Tidal circularization timescale'
tcirc.latex = '\tau_{\rm circ}'
tcirc.label = 'tcirc'
if nplanets eq 0 then tcirc.derive = 0

safronov = parameter
safronov.description = 'Safronov Number'
safronov.latex = '\Theta'
safronov.label = 'safronov'
if nplanets eq 0 then safronov.derive = 0

fave = parameter
fave.unit = '\fluxcgs'
fave.description = 'Incident Flux'
fave.latex = '\fave'
fave.label = 'fave'
fave.cgs = 1d-9
if nplanets eq 0 then fave.derive = 0

u1 = parameter
u1.description = 'linear limb-darkening coeff'
u1.latex = 'u_{1}'
u1.label = 'u1'
u1.scale = 0.15d0
if nplanets eq 0 then u1.derive = 0 $
else u1.fit = 1

u2 = parameter
u2.description = 'quadratic limb-darkening coeff'
u2.latex = 'u_{2}'
u2.label = 'u2'
u2.scale = 0.15d0
if nplanets eq 0 then u2.derive = 0 $
else u2.fit = 1

u3 = parameter
u3.description = 'non-linear limb-darkening coeff'
u3.latex = 'u_{3}'
u3.label = 'u3'
u3.scale = 0.15d0
u3.derive = 0

u4 = parameter
u4.description = 'non-linear limb-darkening coeff'
u4.latex = 'u_{4}'
u4.label = 'u4'
u4.scale = 0.15d0
u4.derive = 0

thermal = parameter
thermal.description = 'Thermal emission from the planet'
thermal.latex = 'A_T'
thermal.label = 'thermal'
thermal.scale = 1d4
thermal.unit = 'ppm'
thermal.fit = 0
thermal.derive = 0

reflect = parameter
reflect.description = 'Reflection from the planet'
reflect.latex = 'A_R'
reflect.label = 'reflect'
reflect.scale = 1d4
reflect.unit = 'ppm'
reflect.fit = 0
reflect.derive = 0

phaseshift = parameter
phaseshift.description = 'Phase shift of reflection curve'
phaseshift.latex = '\phi'
phaseshift.label = 'phaseshift'
phaseshift.scale = 90d0
phaseshift.unit = 'deg'
phaseshift.fit = 0
phaseshift.derive = 0

ellipsoidal = parameter
ellipsoidal.description = 'Ellipsoidal variation of the star'
ellipsoidal.latex = 'A_E'
ellipsoidal.label = 'ellipsoidal'
ellipsoidal.scale = 1d4
ellipsoidal.unit = 'ppm'
ellipsoidal.fit = 0
ellipsoidal.derive = 0

dilute = parameter
dilute.description = 'Dilution from neighboring stars'
dilute.latex = 'A_D'
dilute.label = 'dilute'
dilute.scale = 1d-2
dilute.fit = 0
dilute.derive = 0

mag = parameter
mag.description = 'Apparent Magnitude'
mag.latex = 'mag'
mag.label = 'mag'
mag.scale = 0.1
mag.fit = 0
mag.derive = 0

ttv = parameter
ttv.description = 'Transit Timing Variation'
ttv.latex = 'TTV'
ttv.label = 'ttv'
ttv.unit = 'days'
ttv.scale = 0.02 ;; ~30 minutes
junk = where(fittran, nfittran)
ttv.derive=0

tiv = parameter
tiv.description = 'Transit Inclination Variation'
tiv.latex = 'TiV'
tiv.label = 'tiv'
tiv.unit = 'Radians'
tiv.scale = 0.01*!dpi/180d0
tiv.derive=0

tdeltav = parameter
tdeltav.description = 'Transit Depth Variation'
tdeltav.latex = 'T\delta V'
tdeltav.label = 'tdeltav'
tdeltav.scale = 0.01
tdeltav.derive=0

jitter = parameter
jitter.description = 'RV Jitter'
jitter.latex = '\sigma_J'
jitter.label = 'jitter'
jitter.unit = 'm/s'
jitter.value = 0d0
jitter.scale = 1d0
jitter.fit = 0
if nplanets eq 0 then jitter.derive = 0

jittervar = parameter
jittervar.description = 'RV Jitter Variance'
jittervar.latex = '\sigma_J^2'
jittervar.label = 'jittervar'
jittervar.value = 0d0
jittervar.scale = 1d0
if nplanets eq 0 then jitter.derive = 0 $
else jittervar.fit=1

variance = parameter
variance.description = 'Added Variance'
variance.latex = '\sigma^{2}'
variance.label = 'variance'
variance.value = 0d0
variance.scale = 1d0
if nplanets eq 0 then variance.derive = 0 $
else variance.fit=1

sigma_r = parameter
sigma_r.description = 'Red noise'
sigma_r.latex = '\sigma_r'
sigma_r.label = 'sigma_r'
sigma_r.value = 0d0
sigma_r.scale = 1d0
sigma_r.derive = 0
sigma_r.fit=0

errscale = parameter
errscale.description = 'SED photometry error scaling'
errscale.latex = '\sigma_{SED}'
errscale.label = 'errscale'
errscale.value = 1d0
errscale.scale = 10d0
errscale.fit = 0
errscale.derive = 0

;; Create the structures -- The order here dictates the order in the
;;                          output table.

;; for each star
columnlabels = ''
star = create_struct(mstar.label,mstar,$
                     rstar.label,rstar,$
                     rstarsed.label,rstarsed,$
                     lstar.label,lstar,$
                     fbol.label,fbol,$
                     rhostar.label,rhostar,$
                     logg.label,logg,$
                     teff.label,teff,$
                     teffsed.label,teffsed,$
                     feh.label,feh,$
                     initfeh.label,initfeh,$
                     age.label,age,$
                     eep.label,eep,$
                     logmstar.label,logmstar,$
                     vsini.label,vsini,$
                     vline.label,vline,$                    
                     Av.label,Av,$
                     alpha.label,alpha,$
;                     Ma.label,Ma,$
;                     Mv.label,Mv,$
                     errscale.label,errscale,$
                     ra.label,ra,$       ;; astrometry
                     dec.label,dec,$     ;; astrometry
                     pmra.label,pmra,$   ;; astrometry
                     pmdec.label,pmdec,$ ;; astrometry
                     rvabs.label,rvabs,$ ;; astrometry
                     parallax.label,parallax,$
                     distance.label,distance,$
                     slope.label,slope,$
                     quad.label,quad,$
                     'fluxfile',' ',$
                     'mistsedfile',' ',$
                     'rootlabel','Stellar Parameters:',$
                     'columnlabels',columnlabels,$
                     'label','')
       
;; Unit constants (converted to cgs)
;; As defined by IAU resolutions B2, B3
;; https://www.iau.org/static/resolutions/IAU2012_English.pdf
;; https://arxiv.org/abs/1510.07674, 
;; https://arxiv.org/abs/1507.07956, Table 1
constants = mkconstants()

ndata = 0L
;; if we're fitting an SED, fit the distance, extinction, and error scale
if n_elements(fluxfile) ne 0 then begin
   printandlog, 'WARNING: FLUXFILE has been deprecated. MISTSEDFILE should be used instead.', logname
   printandlog, 'NOTE: When using MISTSEDFILE, the atmosphere is not computed directly or plotted.', logname
   if file_test(fluxfile) then begin
      readcol, fluxfile, junk, format='a', comment='#', /silent
      ndata += n_elements(junk)

      star.fluxfile = fluxfile
      star.errscale.fit = 1
      star.errscale.derive = 1
      star.distance.derive = 1
      star.distance.fit = 1
      star.fbol.derive = 1
      star.parallax.derive = 1
      star.av.fit = 1
      star.av.derive = 1

      if teffsedfloor ne 0d0 then begin
         star.teffsed.fit = 1
         star.teffsed.derive = 1
      endif
      if fbolsedfloor ne 0d0 then begin
         star.rstarsed.fit = 1
         star.rstarsed.derive = 1
      endif

      ;; overwrite the common block, in case it's been called
      ;; before then updated
      ;; the chi2 doesn't matter here, use solar values
      sedchi2 = exofast_sed(fluxfile, 6000d0,1d0,0d0,10d0,logg=4.41d0,met=0d0,alpha=0d0,/redo)
   endif else begin
      printandlog, 'Could not find ' + fluxfile, logname
      stop
   endelse
endif

if n_elements(mistsedfile) ne 0 then begin
   if file_test(mistsedfile) then begin
      readcol, mistsedfile, junk, format='a', comment='#', /silent
      ndata += n_elements(junk)
      star.mistsedfile = mistsedfile
      star.errscale.fit = 1
      star.errscale.derive = 1
      star.distance.derive = 1
      star.distance.fit = 1
      star.fbol.derive = 1
      star.parallax.derive = 1
      star.av.fit = 1
      star.av.derive = 1

      if teffsedfloor ne 0d0 then begin
         star.teffsed.fit = 1
         star.teffsed.derive = 1
      endif
      if fbolsedfloor ne 0d0 then begin
         star.rstarsed.fit = 1
         star.rstarsed.derive = 1
      endif

      ;; overwrite the common block, in case it's been called
      ;; before then updated (without exiting IDL)
      ;; the chi2 doesn't matter here, use solar values
      sedchi2 = mistsed(6000d0, 4.41d0, 0d0, 0d0, 10d0, 1d0, 1d0, mistsedfile, /redo)
   endif else begin
      printandlog, 'Could not find ' + mistsedfile, logname
      stop
   endelse
endif

planet = create_struct($
         period.label,period,$    ;; fundamental (most interesting) parameters
         rp.label,rp,$
         rpearth.label,rpearth,$
         mp.label,mp,$
         mpsun.label,mpsun,$
         logmp.label,logmp,$
         mpearth.label,mpearth,$
         tc.label,tc,$
         tt.label,tt,$
         t0.label,t0,$
         a.label,a,$              
         i.label,i,$
         ideg.label,ideg,$
         e.label,e,$
         omega.label,omega,$
         omegadeg.label,omegadeg,$
         lsinw.label,lsinw,$
         lsinw2.label,lsinw2,$
         lcosw.label,lcosw,$
         bigomega.label,bigomega,$ ;; for astrometry
         bigomegadeg.label,bigomegadeg,$
         lsinbigomega.label,lsinbigomega,$
         lcosbigomega.label,lcosbigomega,$
         teq.label,teq,$
         tcirc.label, tcirc,$
         K.label,k,$              ;; RV parameters
         logK.label,logk,$              ;; RV parameters
         p.label,p,$              ;; Primary Transit parameters
         ar.label,ar,$
         delta.label,delta,$
         depth.label,depth,$
         tau.label,tau,$
         t14.label,t14,$
         tfwhm.label,tfwhm,$
         b.label,b,$
         cosi.label,cosi,$
         bs.label,bs,$          ;; secondary eclipse parameters
         taus.label,taus,$
         t14s.label,t14s,$
         tfwhms.label,tfwhms,$
;         eclipsedepth36.label,eclipsedepth36,$
;         eclipsedepth45.label,eclipsedepth45,$
         eclipsedepth25.label,eclipsedepth25,$
         eclipsedepth50.label,eclipsedepth50,$
         eclipsedepth75.label,eclipsedepth75,$
         rhop.label,rhop,$      ;; less useful parameters            
         rpsun.label,rpsun,$
         logP.label,logp,$  
         loggp.label,loggp,$
         lambda.label,lambda,$
         lambdadeg.label,lambdadeg,$
         lsinlambda.label,lsinlambda,$
         lcoslambda.label,lcoslambda,$
         safronov.label,safronov,$
         fave.label,fave,$
         tp.label,tp,$
         ts.label,ts,$
         ta.label,ta,$
         td.label,td,$
         phase.label,phase,$
         vcve.label,vcve,$
         chord.label,chord,$
         sign.label,sign,$
         sign2.label,sign2,$
         ecosw.label,ecosw,$
         esinw.label,esinw,$
         secosw.label,secosw,$
         sesinw.label,sesinw,$
         qecosw.label,qecosw,$
         qesinw.label,qesinw,$
         msini.label,msini,$
         msiniearth.label,msiniearth,$
         q.label,q,$
         arsun.label,arsun,$
         dr.label,dr,$
         pt.label,pt,$
         ptg.label,ptg,$
         ps.label,ps,$                 
         psg.label,psg,$     
         beam.label,beam,$     ;; other
         'fittran',fittran[0],$        ;; booleans
         'fitrv',fitrv[0],$
         'chen',chen[0],$
         'i180',i180[0],$
         'rossiter',rossiter[0],$
         'fitdt',fitdt[0],$
         'rootlabel','Planetary Parameters:',$
         'label','')

;; for each wavelength
band = create_struct(u1.label,u1,$ ;; linear limb darkening
                     u2.label,u2,$ ;; quadratic limb darkening
                     u3.label,u3,$ ;; 1st non-linear limb darkening
                     u4.label,u4,$ ;; 2nd non-linear limb darkening
                     thermal.label,thermal,$ ;; thermal emission
                     dilute.label,dilute,$   ;; dilution
                     reflect.label,reflect,$ ;; reflection
                     phaseshift.label,phaseshift,$ ;; reflection
                     ellipsoidal.label,ellipsoidal,$
                     eclipsedepth.label,eclipsedepth,$
                     mag.label,mag,$
                     phottobary.label,phottobary,$
                     'name','',$
                     'rootlabel','Wavelength Parameters:',$
                     'label','')

;; for each telescope
telescope = create_struct(gamma.label,gamma,$
                          jitter.label,jitter,$
                          jittervar.label,jittervar,$
                          'rvptrs', ptr_new(),$
                          'detrend',ptr_new(/allocate_heap),$ ;; array of detrending parameters
                          'name','',$
                          'chi2',0L,$
                          'rootlabel','Telescope Parameters:',$
                          'label','')

if ntel le 0 then begin
   telescope.jittervar.fit = 0
   telescope.jittervar.derive = 0
   telescope.jitter.derive = 0
endif

;; for each transit
transit = create_struct(variance.label,variance,$ ;; jitter
                        sigma_r.label,sigma_r,$ ;; Red noise
                        ttv.label,ttv,$           ;; Transit Timing Variation
                        tiv.label,tiv,$ ;; Transit inclination variation
                        tdeltav.label,tdeltav,$ ;; Transit depth variation
                        f0.label,f0,$ ;; normalization
                        'transitptrs',ptr_new(),$ ;; Data
                        'detrend',ptr_new(/allocate_heap),$ ;; array of detrending parameters
                        'bandndx',0L,$
                        'exptime',0d0,$
                        'ninterp',1d0,$
                        'name','',$
                        'epoch',dblarr(nplanets>1) + !values.d_nan,$
                        'pndx',0L,$ ;; index to which planet this corresponds to (-1=>all)
                        'bitmask',0ULL,$ ;; a bitmask for which planet(s) this transit corresponds to
                        'chi2',0L,$
                        'rootlabel','Transit Parameters:',$
                        'rejectflatmodel',0B,$
                        'fitspline',0B,$
                        'splinespace',0.75d0,$
                        'label','') 

doptom = create_struct('dtptrs',ptr_new(),$
                       'rootlabel','Doppler Tomography Parameters:',$
                       'label','',$
                       dtscale.label,dtscale) ;,$
;                       'lambdarange',0,$
;                       'tel','',$
;                       'night','',$
;                       'planetndx')

astrom = create_struct('astromptrs',ptr_new(),$
                       astromscale.label,astromscale,$
                       raoffset.label,raoffset,$
                       decoffset.label,decoffset,$
                       'bandndx',0L,$
                       'rootlabel','Astrometry Parameters:',$
                       'label','',$
                       'epoch',0L)

;; a stellar system has a star, planets, observed bands, observed
;; transits, priors, and global options
ss = create_struct('star',star,$
                   'planet',replicate(planet,nplanets > 1),$
                   'band',replicate(band,nband > 1),$
                   'telescope',replicate(telescope,ntel > 1),$
                   'transit',replicate(transit,ntran > 1),$
                   'doptom',replicate(doptom,ndt>1),$
                   'astrom',replicate(astrom,nastrom>1),$
                   'epochs',dblarr(ntran>1,nplanets>1),$
                   'constants',constants,$
                   'tofit',ptr_new(1),$
                   'priors',ptr_new(1),$
                   'debug',keyword_set(debug),$
                   'verbose',keyword_set(verbose),$
                   'tides',keyword_set(tides),$
                   'randomsign',keyword_set(randomsign),$
                   'ntel',ntel,$
                   'rvepoch',0d0,$
                   'ntran',ntran,$
                   'nastrom',nastrom,$
                   'nband',nband,$
                   'ndt',ndt,$
                   'nplanets',nplanets,$
                   'ndata',ndata,$
                   'mist', ~keyword_set(nomist),$
                   'parsec', keyword_set(parsec),$
                   'fbolsedfloor', fbolsedfloor,$
                   'teffsedfloor', teffsedfloor,$
                   'oned', keyword_set(oned),$
                   'yy', keyword_set(yy),$
                   'torres', keyword_set(torres),$
                   'claret', ~keyword_set(noclaret),$
                   'ttvs', ttvs,$
                   'tivs', tivs,$
                   'tdeltavs', tdeltavs,$
                   'alloworbitcrossing', keyword_set(alloworbitcrossing),$
                   'nsteps',nsteps,$                   
                   'npars',0L,$
                   'burnndx',0L,$
                   'nchains',1L,$
                   'goodchains',1L,$
                   'amoeba',0L,$
                   'logname','',$
                   'chi2',ptr_new(1),$
                   ;; metadata to be able to restart fit
                   'circular', circular,$
                   'fitrv',fitrv,$
                   'fittran',fittran,$
                   'longcadence',longcadence,$
                   'tranpath',tranpath,$
                   'rvpath',rvpath,$
                   'chen',chen,$
                   'planetorder',lindgen(nplanets > 1),$
;                   'dtpath',dtpath,$
;                   'fluxfile',fluxfile,$
                   'delay',delay,$
                   'earth',keyword_set(earth));,$
;                   'prefix',prefix $
;)

if keyword_set(ss.mist) + keyword_set(ss.torres) + keyword_set(ss.yy) gt 1 then begin
   if ~keyword_set(silent) then printandlog, 'WARNING: More than one stellar model invoked -- while this is not forbidden, it is likely an error that will result in underestimated stellar parameters. When using /YY or /TORRES, specify /NOMIST', logname
endif

if n_elements(logname) eq 1 then ss.logname=logname

if ndt gt 0 then begin
   ss.doptom[*].dtptrs = ptrarr(ndt,/allocate_heap)
   if ~keyword_set(silent) then printandlog, 'The index for each DT file is',logname
   for i=0, ndt-1 do begin
      *(ss.doptom[i].dtptrs) = exofast_readdt(dtfiles[i])
      if ~keyword_set(silent) then printandlog, string(i,dtfiles[i],format='(i2,x,a)'),logname
   endfor
   if ~keyword_set(silent) then printandlog, '', logname
endif

;; for each planet 
;; populate the planet fitting parameters
;; planetary labels, a bit optimistic...
plabels = ['b','c','d','e','f','g','h','i','j','k','l','m','n',$
           'o','p','q','r','s','t','u','v','w','x','y','z'] 
for i=0, nplanets-1 do begin
   ss.planet[i].label = plabels[i]

   ;; circular orbit, don't fit e or omega
   if circular[i] then begin
      ss.planet[i].qesinw.fit = 0
      ss.planet[i].qecosw.fit = 0
      ss.planet[i].sesinw.fit = 0
      ss.planet[i].secosw.fit = 0
      ss.planet[i].e.derive=0
      ss.planet[i].omegadeg.derive=0
      ss.planet[i].esinw.derive=0
      ss.planet[i].ecosw.derive=0
   
      ;; don't derive secondary parameters (same as primary if circular)
      ss.planet[i].taus.derive = 0
      ss.planet[i].t14s.derive = 0
      ss.planet[i].tfwhms.derive = 0
      ss.planet[i].bs.derive = 0
      ss.planet[i].ps.derive = 0
      ss.planet[i].psg.derive = 0

   endif
   
   ;; if rossiter is done, fit lambda (for each planet) and vsini (for the star)
   ;; if DT is done, fit lambda (for each planet) and vsini and vline (for the star)
   if rossiter[i] or fitdt[i] then begin
      ss.star.vsini.fit = 1
      ss.star.vsini.derive = 1
      ss.planet[i].lsinlambda.fit = 1
      ss.planet[i].lcoslambda.fit = 1
      ss.planet[i].lambdadeg.derive = 1
      if fitdt[i] then begin
         ss.planet[i].fitdt = 1B
         ss.star.vline.fit = 1
         ss.star.vline.derive = 1
         ss.doptom[*].dtscale.fit = 1
         ss.doptom[*].dtscale.derive = 1
      endif
      if rossiter[i] then ss.planet[i].rossiter = 1B
   endif

   ss.planet[i].fittran = fittran[i]
   ss.planet[i].fitrv = fitrv[i]
   ss.planet[i].chen = chen[i]
   if ss.planet[i].fittran then ss.planet[i].tt.derive = 1B $
   else ss.planet[i].tt.derive = 0B

   if ss.planet[i].fittran and fittt[i] then begin
      ss.planet[i].tt.fit = 1B
      ss.planet[i].tc.fit = 0B
   endif

   ;; fit in Vc/Ve, lcosw, lsinw, and sign if only fitting transit (far more efficient)
   ;; fit in tfwhm0, lcosw, lsinw, and sign if only fitting transit (far more efficient)
   if ~ss.planet[i].fitrv and ss.planet[i].fittran and ~circular[i] and ~novcve[i] then begin
      ss.planet[i].secosw.fit = 0
      ss.planet[i].sesinw.fit = 0
      ss.planet[i].lsinw.fit = 0
      ss.planet[i].lsinw2.fit = 0
      ss.planet[i].lcosw.fit = 0

      ;; reparameterize secosw, sesinw => vcve, lsinw, lcosw
      ss.planet[i].lsinw.fit = 0
      ss.planet[i].lsinw2.fit = 1
      ss.planet[i].lcosw.fit = 1
      ss.planet[i].vcve.fit = 1
      ss.planet[i].vcve.derive = 1
      
      ;; use sign to choose between + and - solutions to the quadratic
      ;; to solve for e (otherwise, use lsinw^2 + lcosw^2 = L)
      if fitsign[i] then begin
         ss.planet[i].sign.fit = 1
         ss.planet[i].sign.derive = 1
;         ss.planet[i].sign2.fit = 1
;         ss.planet[i].sign2.derive = 1
      endif

      ;; reparameterize cosi as chord
      if ~nochord[i] then begin
         ss.planet[i].cosi.fit = 0
         ss.planet[i].chord.fit = 1         
         ss.planet[i].chord.derive = 1
      endif
   endif

   ;; we can marginalize over these parameters 
   ;; even if a transit is not fit.
   ;; And with the Chen & Kipping relation, we can get a decent value
   ;; on the radius and density
   if not fittran[i] and chen[i] then ss.planet[i].cosi.scale = 1d0

   ;; Chen & Kipping prior can be used to fit planetary radius 
   ;; and we can marginalize over cosi
   if not chen[i] then begin
      if not fittran[i] and nastrom eq 0 then begin
         ss.planet[i].cosi.fit = 0
         ss.planet[i].p.fit = 0
         ss.planet[i].cosi.derive = 0
        
         ss.planet[i].p.derive = 0
         ss.planet[i].ideg.derive = 0
         ss.planet[i].delta.derive = 0
         ss.planet[i].depth.derive = 0
         ss.planet[i].mp.derive = 0
         ss.planet[i].mpearth.derive = 0
         
         ;; transit derived pars
         ss.planet[i].b.derive = 0
         ss.planet[i].tfwhm.derive = 0
         ss.planet[i].t14.derive = 0
         ss.planet[i].tau.derive = 0
         ss.planet[i].bs.derive = 0
         
         ;; eclipse derived pars
         ss.planet[i].tfwhms.derive = 0
         ss.planet[i].t14s.derive = 0
         ss.planet[i].taus.derive = 0
         
         ;; requires mass and radius
         ss.planet[i].rhop.derive = 0
         ss.planet[i].loggp.derive = 0
         ss.planet[i].safronov.derive = 0
      endif

   endif

   if fitbeam[i] then begin
      ss.planet[i].beam.fit = 1B
      ss.planet[i].beam.derive = 1B
   endif else if derivebeam[i] then ss.planet[i].beam.derive = 1B

   if i180[i] or ss.nastrom gt 0 then ss.planet[i].i180 = 1

   if fitlogmp[i] then begin
      ss.planet[i].logmp.fit = 1
      ss.planet[i].mpsun.fit = 0
   endif

endfor

for i=0, nband-1 do begin
   ss.band[i].name = bands[i]
   ss.band[i].label = prettybands[(where(bands[i] eq allowedbands))[0]]

   ldcoeffs = quadld(ss.star.logg.value, ss.star.teff.value, ss.star.feh.value, bands[i])
   if finite(ldcoeffs[0]) then ss.band[i].u1.value = ldcoeffs[0] $
   else ss.band[i].u1.value = 0d0
   if finite(ldcoeffs[1]) then ss.band[i].u2.value = ldcoeffs[1] $
   else ss.band[i].u2.value = 0d0

   match = where(fitthermal eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].thermal.fit = 1B
      ss.band[i].thermal.derive = 1B
      ss.band[i].eclipsedepth.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting thermal emission for " + ss.band[i].name + " band",logname
   endif

   match = where(fitreflect eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].reflect.fit = 1B
      ss.band[i].reflect.derive = 1B
      ss.band[i].eclipsedepth.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting reflected light for " + ss.band[i].name + " band", logname
   endif

;   match = where(fitphase eq ss.band[i].name)
;   if match[0] ne -1 then begin
;      ss.band[i].phaseshift.fit = 1B
;      ss.band[i].phaseshift.derive = 1B
;      if ~keyword_set(silent) then printandlog, "Fitting phase offset for " + ss.band[i].name + " band", logname
;   endif

   match = where(fitdilute eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].dilute.fit = 1B
      ss.band[i].dilute.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting dilution for " + ss.band[i].name + " band", logname
   endif

   match = where(fitellip eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].ellipsoidal.fit = 1B
      ss.band[i].ellipsoidal.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting ellipsoidal for " + ss.band[i].name + " band", logname
   endif

endfor

if nband eq 0 then begin
   ss.band[0].u1.fit=0B
   ss.band[0].u2.fit=0B
   ss.band[0].u1.derive=0B
   ss.band[0].u2.derive=0B
endif

;; read in the transit files
if ntran gt 0 then begin
   if n_elements(ninterp) eq 1 and n_elements(exptime) eq 1 then begin
      exptime = dblarr(ntran) + exptime
      ninterp = dblarr(ntran) + ninterp
   endif else if n_elements(ninterp) eq ntran and n_elements(exptime) eq ntran then begin
      ;; do nothing (use exptime and ninterp as is)
   endif else if n_elements(ninterp) eq 0 and n_elements(exptime) eq 0 then begin
      if n_elements(longcadence) eq 0 then begin
         exptime = dblarr(ntran) + 1
         ninterp = dblarr(ntran) + 1
      endif else if n_elements(longcadence) eq 1 then begin
         if longcadence[0] eq 1 then begin
            exptime = dblarr(ntran) + 29.425d0
            ninterp = dblarr(ntran) + 10
         endif else begin
            exptime = dblarr(ntran) + 1
            ninterp = dblarr(ntran) + 1
         endelse
      endif else if n_elements(longcadence) ne ntran then begin
         printandlog, 'LONGCADENCE must be byte or an NTRANSITS (' + strtrim(ntran,2) + ') byte array', logname
         stop
      endif else begin
         exptime = dblarr(ntran) + 29.425d0
         ninterp = dblarr(ntran) + 1
         match = where(longcadence)
         if match[0] ne -1 then ninterp[match] = 10
      endelse
   endif else begin
      printandlog, 'NINTERP and EXPTIME must be unspecified or an NTRANSITS (' + strtrim(ntran,2) + ') array', logname
      stop
   endelse

   if ~keyword_set(silent) then printandlog, 'The index for each transit is',logname
   ss.transit[*].transitptrs = ptrarr(ntran,/allocate_heap)
   for i=0, ntran-1 do begin
      if ~keyword_set(silent) then printandlog, string(i,tranfiles[i],format='(i2,x,a)'),logname
      *(ss.transit[i].transitptrs) = readtran(tranfiles[i], detrendpar=detrend, nplanets=nplanets)

;      ;; create an array of detrending variables 
;      ;; (one for each extra column in the transit file)
      nadd = (*(ss.transit[i].transitptrs)).nadd
;      if nadd ge 1 then *(ss.transit[i].detrendadd) = replicate(detrendadd,nadd)
      nmult = (*(ss.transit[i].transitptrs)).nmult
;      if nmult ge 1 then *(ss.transit[i].detrendmult) = replicate(detrendmult,nmult)

      ss.ndata += n_elements((*(ss.transit[i].transitptrs)).bjd)*(1L+nadd+nmult)

      ss.transit[i].exptime = exptime[i]
      ss.transit[i].ninterp = ninterp[i]
      
      band = (*(ss.transit[i].transitptrs)).band
      ss.transit[i].bandndx = where(ss.band[*].name eq band)
      ss.transit[i].label = (*(ss.transit[i].transitptrs)).label
      ss.transit[i].rejectflatmodel = rejectflatmodel[i]

      ss.transit[i].fitspline = fitspline[i]
      ss.transit[i].splinespace = splinespace[i]
;      ;; F0 and the spline are totally degenerate. 
;      ;; Don't fit (or derive) F0 if flattening with a spline
;      if ss.transit[i].fitspline then begin
;         ss.transit[i].f0.fit = 0B
;         ss.transit[i].f0.derive = 0B
;      endif         

      ;; added variance and added red noise are highly
      ;; degenerate; don't fit both!
      if fitwavelet[i] then begin
         ss.transit[i].variance.fit = 0B
         ss.transit[i].variance.derive = 0B
         ss.transit[i].sigma_r.fit = 1B
         ss.transit[i].sigma_r.derive = 1B
      endif

      for j=0L, ss.nplanets-1L do begin
         if keyword_set(ttvs[i,j]) then begin
            ss.transit[i].ttv.fit = 1B
            ss.transit[i].ttv.derive = 1B
         endif
         if keyword_set(tivs[i,j]) then begin
            ss.transit[i].tiv.fit = 1B
            ss.transit[i].tiv.derive = 1B
         endif
         if keyword_set(tdeltavs[i,j]) then begin
            ss.transit[i].tdeltav.fit = 1B
            ss.transit[i].tdeltav.derive = 1B
         endif
      endfor

   endfor
   if ~keyword_set(silent) then printandlog, '',logname
endif else begin
   ss.transit[0].f0.fit = 0
   ss.transit[0].f0.derive = 0
   ss.transit[0].variance.fit = 0
   ss.transit[0].variance.derive = 0
endelse

;; read in the RV files
if ntel gt 0 then begin
   ss.telescope[*].rvptrs = ptrarr(ntel,/allocate_heap)
   if ~keyword_set(silent) then printandlog, "The index for each RV data set is",logname
   maxpoints = 0
   detrend.scale = 1d0
   for i=0, ntel-1 do begin
      if ~keyword_set(silent) then printandlog, string(i,rvfiles[i],format='(i2,x,a)'),logname
      *(ss.telescope[i].rvptrs) = readrv_detrend(rvfiles[i], detrendpar=detrend)
      ss.telescope[i].label = (*(ss.telescope[i].rvptrs)).label

      ss.ndata += n_elements((*(ss.telescope[i].rvptrs)).bjd)

      ;; create an array of detrending variables 
      ;; (one for each extra column in the rv file)
      nadd = (*(ss.telescope[i].rvptrs)).nadd
      nmult = (*(ss.telescope[i].rvptrs)).nmult
      ss.ndata += n_elements((*(ss.telescope[i].rvptrs)).bjd)*(1L+nadd+nmult)

   endfor
   if ~keyword_set(silent) then printandlog, '', logname
endif else begin
   ss.telescope[*].rvptrs = ptr_new(/allocate_heap)
endelse

for i=0, ss.ntel-1 do begin
   rv = *(ss.telescope[i].rvptrs)
   ss.telescope[i].gamma.value = mean(rv.rv)

;   print, (*ss.telescope[i].rvptrs).planet, (*ss.telescope[i].rvptrs).planet eq -1
;stop

   if (*ss.telescope[i].rvptrs).planet eq -1 then begin
      if n_elements(alltime) eq 0 then begin
         alltime = rv.bjd
         allrv = rv.rv-ss.telescope[i].gamma.value
      endif else begin
         alltime = [alltime,rv.bjd]
         allrv = [allrv,rv.rv-ss.telescope[i].gamma.value]
      endelse
   endif
endfor

;; take a rough stab at gamma, slope, quad, and K
if ss.ntel gt 0 and n_elements(rvepoch) eq 0 then ss.rvepoch = (min(alltime) + max(alltime))/2d0
if ss.star.quad.fit then begin
   coeffs = poly_fit(alltime-ss.rvepoch, allrv, 2, yfit=yfit)
   allrv -= yfit
   ss.star.quad.value = coeffs[2]
   ss.star.slope.value = coeffs[1]
   ss.star.quad.userchanged = 1B
   ss.star.slope.userchanged = 1B
   ss.telescope[*].gamma.value += coeffs[0]
endif else if ss.star.slope.fit then begin
   coeffs = poly_fit(alltime-ss.rvepoch, allrv, 1, yfit=yfit)
   allrv -= yfit
   ss.star.slope.value = coeffs[1]
   ss.telescope[*].gamma.value += coeffs[0]
   ss.star.slope.userchanged = 1B
endif
ss.telescope[*].gamma.userchanged = 1B

if ss.ntel gt 0 then begin
   ss.planet[*].k.value = sqrt(2d0)*stddev(allrv)
   ss.planet[*].k.userchanged=1B
endif

;; read in astrometry files
if nastrom gt 0 then begin
   ss.astrom[*].astromptrs = ptrarr(nastrom,/allocate_heap)
   for i=0L, nastrom-1 do begin
      *(ss.astrom[i].astromptrs) = readastrom(astromfiles[i])
      band = (*(ss.astrom[i].astromptrs)).band
      ss.astrom[i].bandndx = where(ss.band[*].name eq band)
      ss.ndata += n_elements((*(ss.astrom[i].astromptrs)).bjd)
   endfor
endif

for i=0L, nastrom-1 do begin

   astrom = *(ss.astrom[i].astromptrs)
   if ~keyword_set(astrom.userhopa) then begin
      ss.star.parallax.fit = 1B
      ss.star.parallax.derive=1B
      ss.star.ra.fit = 1B
      ss.star.ra.derive=1B
      ss.star.dec.fit = 1B
      ss.star.dec.derive=1B
      ss.star.pmra.fit = 1B
      ss.star.pmra.derive=1B
      ss.star.pmdec.fit = 1B
      ss.star.pmdec.derive=1B
      ss.star.rvabs.fit = 1B
      ss.star.rvabs.derive=1B

      ss.astrom[i].raoffset.fit = 1B
      ss.astrom[i].raoffset.derive=1B
      ss.astrom[i].decoffset.fit = 1B
      ss.astrom[i].decoffset.derive=1B
      ;ss.astrom[i].phottobary.derive = 1B
      ss.astrom[i].astromscale.fit = 1B
      ss.astrom[i].astromscale.derive=1B
   endif
endfor

;; make the prior array
priors = [-1,-1,-1,-1,-1]

;; for each input prior
if ~keyword_set(silent) then printandlog, 'These are the priors that will be applied to the fit.', logname
if ~keyword_set(silent) then printandlog, 'Those with "no prior constraint" only affect the starting values of the fit:', logname
if ~keyword_set(silent) then printandlog, '', logname
openr, lun, priorfile, /get_lun
line = ''
lineno=0
while not eof(lun) do begin
   lineno+=1
   readf, lun, line

   ;; skip commented lines
   if strpos(line,'#') eq 0 then continue

   ;; strip everything after comments
   entries = strsplit((strsplit(line,'#',/extract))[0],/extract)

   nentries = n_elements(entries)
   ;; each line must have at least a name and value
   if nentries lt 2 or nentries gt 6 then begin
      if line ne '' and ~keyword_set(silent) then printandlog, 'WARNING: line ' + strtrim(lineno,2) + ' in ' + priorfile + ' is not legal syntax (NAME VALUE [UNCERTAINTY] [LOWERBOUND] [UPPERBOUND] [START]); ignoring: ' + line, logname
      continue
   endif 

   ;; extract Name, value, uncertainty, lowerbound, and upper bound
   ;; (or default to Name, Value, no uncertainty, -Inf, +Inf)
   priorname = entries[0]

   ;; for MIST models, age is derived after the fact, which wouldn't allow a prior. 
   ;; So if a prior is supplied, make sure we fit it.
   if strupcase(priorname) eq 'AGE' and n_elements(entries) gt 2 then ss.star.age.fit = 1B

   priorval = double(entries[1])
   if nentries ge 3 then priorwidth = double(entries[2]) $
   else priorwidth = -1
   if priorwidth gt 0d0 then ss.ndata++
   if nentries ge 4 then begin
      if strupcase(entries[3]) eq '-INF' then lowerbound = -!values.d_infinity $
      else lowerbound = double(entries[3])
   endif else lowerbound = -!values.d_infinity
   if nentries ge 5 then begin
      if strupcase(entries[4]) eq 'INF' then upperbound = !values.d_infinity $
      else upperbound = double(entries[4])
   endif else upperbound = !values.d_infinity
   if nentries ge 6 then begin
      startval = double(entries[5])
   endif else startval = priorval

   ;; determine the subscript
   tmp = strsplit(priorname,'_',/extract)
   if n_elements(tmp) eq 2 then priornum = tmp[1] $
   else priornum = '0'
   priorlabel = tmp[0]

   ;; allow labels 'b','c','d', etc.
   numndx = (where(plabels eq strtrim(priornum,2)))[0]
   if numndx ne -1 then priornum = numndx
   priornum = long(priornum)

   found=0
   ;; look for the name in the structure
   for i=0L, n_tags(ss)-1 do begin
      if (n_elements(ss.(i))-1) ge priornum then begin
         for k=0, n_tags(ss.(i)[priornum])-1 do begin

            ;; this captures the detrending variables
            if (size(ss.(i)[priornum].(k)))[1] eq 10 then begin ;; if it's a pointer
               if ptr_valid(ss.(i)[priornum].(k)) then begin
                  for l=0L, n_tags(*(ss.(i)[priornum].(k)))-1 do begin

                     if (size((*(ss.(i)[priornum].(k))).(l)))[2] eq 8 then begin

                        valid = 0
                        if strpos(strupcase(priorlabel),'C') eq 0 then begin
                           if strmid((*(ss.(i)[priornum].(k))).(l)[0].label,0,1) eq 'C' then begin ;; if it's the additive variable
                              if (*(ss.(i)[priornum].(k))).nadd gt 0 then begin
                                 detrendnum = long((strsplit(strupcase(priorlabel),'C',/extract))[0])
                                 if detrendnum lt (*(ss.(i)[priornum].(k))).nadd then valid = 1
                              endif
                           endif
                        endif else if strpos(strupcase(priorlabel),'M') eq 0 then begin
                           if strmid((*(ss.(i)[priornum].(k))).(l)[0].label,0,1) eq 'M' then begin
                              if (*(ss.(i)[priornum].(k))).nmult gt 0 then begin
                                 detrendnum = long((strsplit(strupcase(priorlabel),'M',/extract))[0])
                                 if detrendnum lt (*(ss.(i)[priornum].(k))).nmult then valid = 1
                              endif
                           endif
                        endif
                           
                        if valid then begin

                           if (*(ss.(i)[priornum].(k))).(l)[detrendnum].label eq strupcase(priorlabel) then begin

                              if not (*(ss.(i)[priornum].(k))).(l)[detrendnum].fit and not (*(ss.(i)[priornum].(k))).(l)[detrendnum].derive then begin
                                 ;; it's not fit or derived; 
                                 ;; only use it to derive the starting parameters
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].value = startval
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].userchanged = 1B

                                 if nplanets ne 0 then $
                                    if ~keyword_set(silent) then printandlog, "WARNING: Prior supplied on '" + $
                                      priorname + "' but it is neither fitted or derived. Not applying prior, " + $
                                      'but it will be used to derive the starting parameters if possible.', logname
                                 
                              endif else begin
                              
                                 ;; found it! change the default starting guess to the value
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].prior = priorval
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].value = startval
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].upperbound = upperbound
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].lowerbound = lowerbound
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].userchanged = 1B
                                 
                                 if priorwidth eq 0d0 then begin
                                    if ~(*(ss.(i)[priornum].(k))).(l)[detrendnum].fit then begin
                                       if ~keyword_set(silent) then printandlog, priorname + ' is not fit. Only fitted parameters can be fixed. Ignoring constraint', logname
                                    endif else begin
                                       ;; priorwidth = 0 => fix it at the prior value
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].fit = 0d0
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].derive = 0d0
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].priorwidth = 0d0

                                       if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (fixed)', logname
                                    endelse
                                 endif else if finite(priorwidth) and priorwidth gt 0d0 or finite(lowerbound) or finite(upperbound) then begin
                                    ;; apply a Gaussian prior with width = priorwidth
                                    priors=[[priors],[i,priornum,k,l,detrendnum]]
                                    
                                    if priorwidth gt 0 and finite(priorwidth) then begin
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].scale = priorwidth*3d0
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].priorwidth = priorwidth
                                    endif
                                    
                                    if ~keyword_set(silent) then begin
                                       if priorwidth lt 0d0 then printandlog, priorname + ' = ' + strtrim(startval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname $
                                       else begin
                                          printandlog, priorname + ' = ' + strtrim(priorval,2) + ' +/- ' + strtrim(priorwidth,2) + '; bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                                          if startval ne priorval then begin
                                             printandlog, 'NOTE: ' + priorname + ' starts at ' + strtrim(startval,2), logname
                                          endif
                                       endelse
                                    endif
                                    
                                 endif else begin
                                    ;; else no prior, just change the default starting value
                                    if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(startval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                                 endelse
                                 
                              endelse
                              found=1
                              break
                           endif ; else found=0
                           
                           if found then break
                           
                        endif
                        
                     endif
                     
                  endfor
               endif          
            endif else if n_tags(ss.(i)[priornum].(k)) ne 0 then begin
               ;; and this captures everything else
               if tag_exist(ss.(i)[priornum],priorlabel,index=ndx) then begin
                  
                  circ = 0
                  if tag_exist(ss.(i),'e') then if circular[priornum] then circ = 1
                  
                  ;; if circular flag set, ignore priors on any eccentricity-derived parameter
                  if ss.(i)[priornum].(ndx).label eq 'e' or ss.(i)[priornum].(ndx).label eq 'omega' or $
                     ss.(i)[priornum].(ndx).label eq 'omegadeg' or $
                     ss.(i)[priornum].(ndx).label eq 'secosw' or ss.(i)[priornum].(ndx).label eq 'sesinw' or $
                     ss.(i)[priornum].(ndx).label eq 'ecosw' or ss.(i)[priornum].(ndx).label eq 'esinw' and circ then begin                        
                     if ~keyword_set(silent) then printandlog, "WARNING: Prior supplied on '" + $
                         priorname + "' but planet is specified to be circular; not applying prior", logname
                  endif else if not ss.(i)[priornum].(ndx).fit and not ss.(i)[priornum].(ndx).derive then begin
                     ;; it's not fit or derived; 
                     ;; only use it to derive the starting parameters
                     ss.(i)[priornum].(ndx).value = startval
                     ss.(i)[priornum].(ndx).userchanged = 1B
                     if nplanets ne 0 then $
                        if ~keyword_set(silent) then printandlog, "WARNING: Prior supplied on '" + $
                            priorname + "' but it is neither fitted or derived. Not applying prior, " + $
                            'but it will be used to derive the starting parameters if possible.', logname
                  endif else begin
                     
                     ;; found it! change the default starting guess to the value
                     ss.(i)[priornum].(ndx).prior = priorval
                     ss.(i)[priornum].(ndx).value = startval
                     ss.(i)[priornum].(ndx).upperbound = upperbound
                     ss.(i)[priornum].(ndx).lowerbound = lowerbound
                     ss.(i)[priornum].(ndx).userchanged = 1B
                     
                     if priorwidth eq 0d0 then begin
                        if ~ss.(i)[priornum].(ndx).fit then begin
                           if ~keyword_set(silent) then printandlog, priorname + ' is not fit. Only fitted parameters can be fixed. Ignoring constraint', logname
                        endif else begin
                           ;; priorwidth = 0 => fix it at the prior value
                           ss.(i)[priornum].(ndx).fit = 0d0
                           ss.(i)[priornum].(ndx).derive = 0d0
                           ss.(i)[priornum].(ndx).priorwidth = 0d0
                           if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (fixed)', logname
                        endelse
                     endif else if finite(priorwidth) and priorwidth gt 0d0 or finite(lowerbound) or finite(upperbound) then begin
                        ;; apply a Gaussian prior with width = priorwidth
                        priors=[[priors],[i,priornum,ndx,-1,-1]]
                        
                        if priorwidth gt 0 and finite(priorwidth) then begin
                           ss.(i)[priornum].(ndx).scale = priorwidth*3d0
                           ss.(i)[priornum].(ndx).priorwidth = priorwidth
                        endif                  
                        
                        if ~keyword_set(silent) then begin
                           if priorwidth lt 0 then printandlog, priorname + ' = ' + strtrim(startval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname $
                           else begin
                              printandlog, priorname + ' = ' + strtrim(priorval,2) + ' +/- ' + strtrim(priorwidth,2) + '; bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                              if startval ne priorval then begin
                                 printandlog, 'NOTE: ' + priorname + ' starts at ' + strtrim(startval,2), logname
                              endif
                           endelse
                        endif
                        
                     endif else begin
                        ;; else no prior, just change the default starting value
                        if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(startval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                     endelse
                  endelse

                  found=1
                  i = n_tags(ss)-1L
                  break
               endif; else found=0
            endif

         endfor
      endif
   endfor

   ;; didn't find it, warn user
   if not found and nplanets ne 0 and ~keyword_set(silent) then printandlog, "WARNING: No parameter matches '" + $
      priorname + "' from " + priorfile + "; not applying prior", logname
                
endwhile

free_lun, lun
if ~keyword_set(silent) then printandlog, '', logname

;; do we have enough information to derive the distance?
;if (where(priorname eq 'distance'))[0] ne -1 

priors = priors[*,1:*]
*(ss.priors) = priors

;; determine the epoch for each observation
for i=0, ntran-1 do begin

;   if i eq 0 then begin
;      ;; variations are defined relative to the first transit
;      if keyword_set(tivs) then begin
;         ss.transit[i].tiv.fit = 0
;         ss.transit[i].tiv.derive = 0
;      endif
;      if keyword_set(tdeltavs) then begin
;         ss.transit[i].tdeltav.fit = 0
;         ss.transit[i].tdeltav.derive = 0
;      endif
;   endif
   
   for j=0, nplanets-1 do begin

      

;      ss.transit[i].epoch = min(round((mean((*(ss.transit[i].transitptrs)).bjd) - ss.planet[ss.transit[i].pndx].tc.value)/ss.planet[ss.transit[i].pndx].period.value))

      if n_elements(nvalues) ne 0 then begin
         tc = median(ss.planet[j].tc.value)
         period = median(ss.planet[j].period.value)
         if period eq 0 then period = median(10^ss.planet[j].logp.value)
      endif else begin
         tc = ss.planet[j].tc.value
         period = ss.planet[j].period.value
         if period eq 0 then period = 10^ss.planet[j].logp.value
      endelse
      
      minbjd = min((*(ss.transit[i].transitptrs)).bjd, max=maxbjd)
      

      ss.epochs[i,j] = round(((maxbjd+minbjd)/2d0 - tc)/period)

      
;      t14 = ??
;      t14s = ??

;      thistc = tc + epoch*period
;      thists = 
;      if (minbjd le thistc + t14 and maxbjd ge thistc - t14) or $
;         (minbjd le thists + t14s and maxbjd ge thists - t14s) then $
;            ;; the ith transit covers the jth planet; add it to the bitmask
;            transit[i].bitmask += 2ULL^j 
;

      ss.transit[i].epoch[j] = ss.epochs[i,j]

   endfor
endfor


;; creates an array of indicies into the stellar structure to map which parameters should be fit
;; fit[*,0] indexes the object [star=0, planet=1, band=2, telescope=3, or transit=4]
;; fit[*,1] each object can have any number of copies. 
;; This indexes which copy (e.g., planet b=0, planet c=1 or B band=0, V band=1) 
;; fit[*,2] indexes the parameter of the object (e.g., Teff=0, [Fe/H]=1)
;; this assumes a certain structure of the parameters... is that ok?
tofit = [-1,-1,-1,-1,-1]
for i=0L, n_tags(ss)-1 do begin
   for j=0L, n_elements(ss.(i))-1 do begin
      for k=0L, n_tags(ss.(i)[j])-1 do begin

         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin
            if ptr_valid(ss.(i)[j].(k)) then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'fit') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].fit then tofit = [[tofit],[i,j,k,l,m]]
                           if (*(ss.(i)[j].(k))).(l)[m].fit or (*(ss.(i)[j].(k))).(l)[m].derive then ss.npars++
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            ;; and this captures everything else
            if tag_exist(ss.(i)[j].(k),'fit') then begin
               if ss.(i)[j].(k).fit then tofit = [[tofit],[i,j,k,-1,-1]]
               if ss.(i)[j].(k).fit or ss.(i)[j].(k).derive then ss.npars++
            endif
         endif
      endfor
   endfor
endfor
tofit = tofit[*,1:*]
*(ss.tofit) = tofit
ss.nchains = n_elements((*ss.tofit)[0,*])*2L

;; populate the best values
if arg_present(best) then pars2str, best, ss, /best

;; don't do these when creating the MCMC structure
if n_elements(ss.star.mstar.value) eq 1 then begin

   ;; derive all step parameters
   if not pars2step(ss) then begin
      printandlog, 'ERROR: The isochrones are not applicable here; refine your priors.', logname
      stop
   endif

   ;; use this to require planets to stay in the same order as they start
   if ss.nplanets gt 0 then begin
      priorsupplied = where(ss.planet.period.value ne 0d0)
      if priorsupplied[0] ne -1 then ss.planet[priorsupplied].logp.value = alog10(ss.planet[priorsupplied].period.value)
      ss.planetorder = sort(ss.planet.logp.value)
   endif

   ;; return an error if the starting stellar system is not allowed
   if step2pars(ss,/verbose,/changedefaults) eq -1 then begin
      printandlog, 'ERROR: starting values for stellar system not allowed; refine priors', logname
      stop
   endif

   ;; calculate the right period stepping scale for AMOEBA based on the
   ;; input data
   perscale, ss
endif else begin
   ;; use this to require planets to stay in the same order as they start
   if ss.nplanets gt 0 then ss.planetorder = sort(ss.planet.period.value[0])
endelse

;; load stellar structure into common block
;; must do it here because threads can't pass structures
if n_elements(chi2func) eq 1 then junk = call_function(chi2func, /loadss, ss0=ss)

if (ss.star.mistsedfile ne ' ' or ss.star.fluxfile ne ' ') and $
   (nastrom eq 0 and ~finite(ss.star.parallax.priorwidth)) then  begin
   printandlog, 'WARNING: Fitting an SED without providing parallax information will essentially determine a photometric parallax.', logname
endif

if keyword_set(ss.mist) then ss.ndata+=3
if keyword_set(ss.yy) then ss.ndata+=3
if keyword_set(ss.parsec) then ss.ndata+=3
if keyword_set(ss.torres) then ss.ndata+=2

return, ss

end
