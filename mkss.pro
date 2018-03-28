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
               fitthermal=fitthermal, fitreflect=fitreflect, fitdilute=fitdilute,$
               nvalues=nvalues, debug=debug, verbose=verbose, priorfile=priorfile, $
               rvpath=rvpath, tranpath=tranpath, dtpath=dtpath, fluxfile=fluxfile, $
               longcadence=longcadence, ninterp=ninterp, exptime=exptime,$
               earth=earth, silent=silent, noyy=noyy, torres=torres, mist=mist, $
               noclaret=noclaret,alloworbitcrossing=alloworbitcrossing,$
               logname=logname, best=best

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

if not keyword_set(longcadence) then longcadence=0B
if n_elements(fitthermal) eq 0 then fitthermal = []
if n_elements(fitreflect) eq 0 then fitreflect = []
if n_elements(fitdilute) eq 0 then fitdilute = []
if n_elements(tranpath) eq 0 then tranpath = ''

;; read in the transit files
if tranpath ne '' then begin
   tranfiles=file_search(tranpath,count=ntran)
   if ntran eq 0 then begin
      printandlog, "No transit files files found matching " + strtrim(tranpath,2) + "; please check TRANPATH", logname
      stop
   endif

   ;; find the unique bands
   bands = tranfiles
   for i=0, ntran-1 do begin
      bands[i] = (strsplit(tranfiles[i],'.',/extract))(1)
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
      allowedbands = ['U','B','V','R','I','J','H','K',$
                      'Sloanu','Sloang','Sloanr','Sloani','Sloanz',$
                      'Kepler','TESS','CoRoT','Spit36','Spit45','Spit58','Spit80',$
                      'u','b','v','y']
      prettybands = ['U','B','V','R','I','J','H','K',$
                     "u'","g'","r'","i'","z'",$
                     'Kepler','TESS','CoRoT','$3.6\mu m$','$4.5\mu m$','$5.8\mu m$','$8.0\mu m$',$
                     'u','b','v','y']
      if (where(bands eq bands[i]))[0] eq -1 then begin
         printandlog, 'ERROR: band (' + bands[i] + ') not recognized; please select one of the following:'
         printandlog, string(allowedbands)
         stop
      endif
   endfor

   bands = bands[uniq(bands, sort(bands))]
   nband = n_elements(bands)
   printandlog, 'The index for each fitted band is', logname
   for i=0L, nband-1 do printandlog, string(i,bands[i],format='(i2,x,a)'),logname
   printandlog, '', logname
endif else begin
   ntran = 0
   nband = 0
   tranpath = ''
endelse

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

if n_elements(circular) ne nplanets or nplanets eq 0 then circular = bytarr(nplanets>1)
if n_elements(rossiter) ne nplanets or nplanets eq 0 then rossiter = bytarr(nplanets>1)
if n_elements(fitdt) ne nplanets or nplanets eq 0 then fitdt = bytarr(nplanets>1)

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
   if (where((~fitrv) and (~fittran)))[0] ne -1 then begin
      printandlog, 'Either a transit or RV must be fit for each planet', logname
      stop
   end
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
                          'best',!values.d_nan,$ ;; best-fit parameter
                          'latex','',$        ;; latex label for the latex table
                          'description','',$  ;; a verbose description
                          'unit','',$         ;; units ("Radians" triggers special angular behavior)
                          'medvalue','',$     ;; median value
                          'upper','',$        ;; upper error bar (68% confidence)
                          'lower','',$        ;; lower error bar (68% confidence)
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
teff.value = 6000d0
teff.unit = 'K'
teff.description = 'Effective Temperature'
teff.latex = 'T_{\rm eff}'
teff.label = 'teff'
teff.fit=1
teff.scale = 500d0

feh = parameter
feh.value = 0d0
feh.description = 'Metalicity'
feh.latex = '[{\rm Fe/H}]'
feh.label = 'feh'
feh.fit=1
feh.scale = 0.5d0

initfeh = parameter
initfeh.value = 0d0
initfeh.description = 'Initial Metalicity'
initfeh.latex = '[{\rm Fe/H}]_{0}'
initfeh.label = 'initfeh'
initfeh.scale = 0.5d0
if keyword_set(mist) then begin
   initfeh.fit=1
   initfeh.derive=1
endif else begin
   initfeh.fit=0
   initfeh.derive=0
endelse

Av = parameter
Av.description = 'V-band extinction'
Av.latex = 'A_v'
Av.label = 'Av'
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
distance.fit = 0
distance.derive = 0
distance.scale = 100
distance.value = 10

parallax = parameter
parallax.unit = 'mas'
parallax.description = 'Parallax'
parallax.latex = '\pi'
parallax.label = 'parallax'
parallax.cgs = 3600d3*180d0/!dpi ;; rad/mas
parallax.derive = 0

gamma = parameter
gamma.unit = 'm/s'
gamma.description = 'Instrumental offset'
gamma.latex = '\gamma'
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
if keyword_set(fitslope) then slope.fit = 1 $
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

age = parameter
age.value = 4.6d0
age.unit = 'Gyr'
age.description = 'Age'
age.latex = 'Age'
age.label = 'age'
age.cgs = 3600d0*24d0*365.242d0*1d9
age.derive = 0
if ~keyword_set(noyy) then begin
   age.fit = 1
   age.derive = 1
endif
if keyword_set(mist) then age.derive=1

age.scale = 3d0

eep = parameter
eep.value = 355.65d0 ;; solar EEP
eep.unit = ''
eep.description = 'Equal Evolutionary Point'
eep.latex = 'EEP'
eep.label = 'eep'
eep.scale = 300d0 
if keyword_set(mist) then begin
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
tc.description = 'Time of Transit'
tc.latex = 'T_C'
tc.label = 'tc'
tc.cgs = 86400d0
if nplanets eq 0 then tc.derive=0 $
else tc.fit = 1
tc.scale = 0.1

t0 = parameter
t0.unit = '\bjdtdb'
t0.description = 'Optimal Transit Time'
t0.latex = 'T_0'
t0.label = 't0'
t0.cgs = 86400d0
if nplanets eq 0 then t0.derive=0

logp = parameter
logp.unit = ''
logp.description = 'Log of Period'
logp.latex = 'log P'
logp.label = 'logp'
logp.cgs = !values.d_nan
if nplanets ne 0 then logp.fit = 1
logp.derive=0
logp.scale = 0.01d0

qesinw = parameter
qesinw.unit = ''
qesinw.description = ''
qesinw.latex = 'e^{1/4} sin{\omega_*}'
qesinw.label = 'qesinw'
qesinw.cgs = !values.d_nan
qesinw.scale = 0.1d0
qesinw.derive = 0

qecosw = parameter
qecosw.unit = ''
qecosw.description = ''
qecosw.latex = 'e^{1/4} cos{\omega_*}'
qecosw.label = 'qecosw'
qecosw.cgs = !values.d_nan
qecosw.scale = 0.1d0
qecosw.derive = 0

sesinw = parameter
sesinw.unit = ''
sesinw.description = ''
sesinw.latex = '\sqrt{e} sin{\omega_*}'
sesinw.label = 'sesinw'
sesinw.cgs = !values.d_nan
sesinw.scale = 0.1d0
sesinw.derive = 0

secosw = parameter
secosw.unit = ''
secosw.description = ''
secosw.latex = '\sqrt{e} cos{\omega_*}'
secosw.label = 'secosw'
secosw.cgs = !values.d_nan
secosw.scale = 0.1d0
secosw.derive = 0

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
esinw.latex = 'esin{\omega_*}'
esinw.label = 'esinw'
esinw.cgs = !values.d_nan
if nplanets eq 0 then esinw.derive=0

ecosw = parameter
ecosw.unit = ''
ecosw.description = ''
ecosw.latex = 'ecos{\omega_*}'
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

omegadeg = parameter
omegadeg.unit = 'Degrees'
omegadeg.description = 'Argument of Periastron'
omegadeg.latex = '\omega_*'
omegadeg.label = 'omegadeg'
omegadeg.cgs = !dpi/180d0
if nplanets eq 0 then omegadeg.derive=0

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
cosi.latex = 'cos i'
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

eclipsedepth36 = parameter
eclipsedepth36.unit = 'ppm'
eclipsedepth36.description = 'Blackbody eclipse depth at 3.6$\mu$m'
eclipsedepth36.latex = '\delta_{S,3.6\mu m}'
eclipsedepth36.label = 'eclipsedepth36'
eclipsedepth36.cgs = 1d6
if nplanets eq 0 then eclipsedepth36.derive = 0

eclipsedepth45 = parameter
eclipsedepth45.unit = 'ppm'
eclipsedepth45.description = 'Blackbody eclipse depth at 4.5$\mu$m'
eclipsedepth45.latex = '\delta_{S,4.5\mu m}'
eclipsedepth45.label = 'eclipsedepth45'
eclipsedepth45.cgs = 1d6
if nplanets eq 0 then eclipsedepth45.derive = 0

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

;sqrtcosicoslambda = parameter
;sqrtcosicoslambda.description = 'sqrtcosicoslambda'
;sqrtcosicoslambda.latex = '\sqrt{cosi_*}\cos{\lambda}'
;sqrtcosicoslambda.label = 'sqrtcosicoslambda'
;sqrtcosicoslambda.scale = 1d2
;sqrtcosicoslambda.derive = 0
;
;sqrtsinisinlambda = parameter
;sqrtsinisinlambda.description = 'sqrtsinisinlambda'
;sqrtsinisinlambda.latex = '\sqrt{sini}\sin{\lambda}'
;sqrtsinisinlambda.label = 'sqrtsinisinlambda'
;sqrtsinisinlambda.scale = 1d2
;sqrtsinisinlambda.derive = 0

lambda = parameter
lambda.unit = 'Radians'
lambda.description = 'Projected Spin-orbit alignment'
lambda.latex = '\lambda'
lambda.label = 'lambda'
lambda.scale = !dpi
lambda.derive = 0

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

macturb = parameter
macturb.unit = 'm/s'
macturb.description = 'Macroturbulence'
macturb.latex = 'macturb'
macturb.label = 'macturb'
macturb.cgs = 1000d0
macturb.derive = 0
macturb.value = 4d3
macturb.scale = 2d3

dtscale = parameter
dtscale.description = 'Doppler Tomography Error scaling'
dtscale.latex = '\sigma_{DT}'
dtscale.label = 'dtscale'
dtscale.derive = 0
dtscale.scale = 1d2
dtscale.value = 1d0

logk = parameter
logk.description = 'Log of RV semi-amplitude'
logk.latex = 'logK'
logk.label = 'logk'
if nplanets eq 0 then logk.derive=0 $
else logk.fit = 1
logk.value = 1d0
logk.scale = 1d0

k = parameter
k.value = 10^logk.value ;; earth
k.unit = 'm/s'
k.description = 'RV semi-amplitude'
k.latex = 'K'
k.label = 'k'
k.cgs = 100d0
if nplanets eq 0 then k.derive=0

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
beam.description = 'Doppler Beaming'
beam.latex = 'A_B'
beam.label = 'beam'
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
if nplanets eq 0 then mpsun.derive = 0

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
if keyword_set(ttvs) then begin
   ttv.fit = 1
endif else ttv.derive=0

tiv = parameter
tiv.description = 'Transit Inclination Variation'
tiv.latex = 'TIV'
tiv.label = 'tiv'
tiv.unit = 'Radians'
tiv.scale = 0.01*!dpi/180d0
if keyword_set(tivs) then begin
   tiv.fit = 1
endif else tiv.derive=0

tdeltav = parameter
tdeltav.description = 'Transit Depth Variation'
tdeltav.latex = 'T\delta V'
tdeltav.label = 'tdeltav'
tdeltav.scale = 0.01
if keyword_set(tdeltavs) then begin
   tdeltav.fit = 1
endif else tdeltav.derive=0

rednoise = parameter
rednoise.description = 'Red Noise'
rednoise.latex = '\sigma_r'
rednoise.label = 'rednoise'
rednoise.fit=0
rednoise.derive=0

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
                     lstar.label,lstar,$
                     rhostar.label,rhostar,$
                     logg.label,logg,$
                     teff.label,teff,$
                     feh.label,feh,$
                     initfeh.label,initfeh,$
                     age.label,age,$
                     eep.label,eep,$
                     logmstar.label,logmstar,$
                     vsini.label,vsini,$
                     macturb.label,macturb,$
                     Av.label,Av,$
                     alpha.label,alpha,$
;                     Ma.label,Ma,$
;                     Mv.label,Mv,$
                     errscale.label,errscale,$
                     distance.label,distance,$
                     parallax.label,parallax,$
;                     ellip.label,0d0,$
;                     ra.label,0d0,$       ;; for astrometry?
;                     dec.label,0d0,$ ;; for astrometry?
;                     bigomega.label,0d0,$ ;; for astrometry
;                     pmra.label,pmra,$    ;; astrometry
;                     pmdec.label,pmdec,$ ;; astrometry
                     slope.label,slope,$
                     quad.label,quad,$
                     'fluxfile','',$
                     'rootlabel','Stellar Parameters:',$
                     'columnlabels',columnlabels,$
                     'label','')
       
;; Unit constants (converted to cgs)
;; As defined by IAU resolutions B2, B3
;; https://www.iau.org/static/resolutions/IAU2012_English.pdf
;; https://arxiv.org/abs/1510.07674, 
;; https://arxiv.org/abs/1507.07956, Table 1
constants = mkconstants()

;; if we're fitting an SED, fit the distance, extinction, and error scale
if n_elements(fluxfile) ne 0 then begin
   if file_test(fluxfile) then begin
      star.fluxfile = fluxfile
      star.errscale.fit = 1
      star.errscale.derive = 1
      star.distance.fit = 1
      star.distance.derive = 1
      star.parallax.derive = 1
      star.av.fit = 1
      star.av.derive = 1
   endif else printandlog, 'Could not find ' + fluxfile, logname
endif

planet = create_struct($
         period.label,period,$    ;; fundamental (most interesting) parameters
         rp.label,rp,$
         rpearth.label,rpearth,$
         tc.label,tc,$
         t0.label,t0,$
         a.label,a,$              
         i.label,i,$
         ideg.label,ideg,$
         e.label,e,$
         omega.label,omega,$
         omegadeg.label,omegadeg,$
         teq.label,teq,$
         mp.label,mp,$
         mpsun.label,mpsun,$
         mpearth.label,mpearth,$
         K.label,k,$              ;; RV parameters
         logK.label,logk,$
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
         eclipsedepth36.label,eclipsedepth36,$
         eclipsedepth45.label,eclipsedepth45,$
         rhop.label,rhop,$      ;; less useful parameters            
         rpsun.label,rpsun,$
         logP.label,logp,$  
         loggp.label,loggp,$
         lambda.label,lambda,$
         lambdadeg.label,lambdadeg,$
         safronov.label,safronov,$
         fave.label,fave,$
         tp.label,tp,$
         ts.label,ts,$
         ta.label,ta,$
         td.label,td,$
         phase.label,phase,$
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
                     mag.label,mag,$
                     'name','',$
                     'rootlabel','Wavelength Parameters:',$
                     'label','')

;; for each telescope
telescope = create_struct(gamma.label,gamma,$
                          jitter.label,jitter,$
                          jittervar.label,jittervar,$
                          'rvptrs', ptr_new(),$
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
transit = create_struct(variance.label,variance,$ ;; Red noise
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
                        'chi2',0L,$
                        'rootlabel','Transit Parameters:',$
                        'label','') 

if n_elements(dtpath) ne 0 then begin
   dtfiles = file_search(dtpath)
endif else dtfiles = []
ndt = n_elements(dtfiles)
doptom = create_struct('dtptrs',ptr_new(),$
                       'rootlabel','Doppler Tomography Parameters:',$
                       'label','',$
                       dtscale.label,dtscale) ;,$
;                       'lambdarange',0,$
;                       'tel','',$
;                       'night','',$
;                       'planetndx')

;; a stellar system has a star, planets, observed bands, observed
;; transits, priors, and global options
ss = create_struct('star',star,$
                   'planet',replicate(planet,nplanets > 1),$
                   'band',replicate(band,nband > 1),$
                   'telescope',replicate(telescope,ntel > 1),$
                   'transit',replicate(transit,ntran > 1),$
                   'doptom',replicate(doptom,ndt>1),$
                   'constants',constants,$
                   'tofit',ptr_new(1),$
                   'priors',ptr_new(1),$
                   'debug',keyword_set(debug),$
                   'verbose',keyword_set(verbose),$
                   'tides',keyword_set(tides),$
                   'ntel',ntel,$
                   'ntran',ntran,$
                   'nband',nband,$
                   'ndt',ndt,$
                   'nplanets',nplanets,$
                   'mist', keyword_set(mist),$
                   'yy', ~keyword_set(noyy),$
                   'torres', keyword_set(torres),$
                   'claret', ~keyword_set(noclaret),$
                   'ttvs', keyword_set(ttvs),$
                   'tivs', keyword_set(tivs),$
                   'tdeltavs', keyword_set(tdeltavs),$
                   'alloworbitcrossing', keyword_set(alloworbitcrossing),$
                   'nsteps',nsteps,$
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
;                   'dtpath',dtpath,$
;                   'fluxfile',fluxfile,$
                   'earth',keyword_set(earth));,$
;                   'prefix',prefix $
;)

if n_elements(logname) eq 1 then ss.logname=logname

if ndt gt 0 then begin
   ss.doptom[*].dtptrs = ptrarr(ndt,/allocate_heap)
   for i=0, n_elements(dtfiles)-1 do begin
      *(ss.doptom[i].dtptrs) = exofast_readdt(dtfiles[i])
   endfor
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
   ;; if DT is done, fit lambda (for each planet) and vsini and macturb (for the star)
   if rossiter[i] or fitdt[i] then begin
      ss.star.vsini.fit = 1
      ss.star.vsini.derive = 1
      ss.planet[i].lambda.fit = 1
      ss.planet[i].lambdadeg.derive = 1
      if fitdt[i] then begin
         ss.planet[i].fitdt = 1B
         ss.star.macturb.fit = 1
         ss.star.macturb.derive = 1
         ss.doptom[i].dtscale.fit = 1
         ss.doptom[i].dtscale.derive = 1
      endif
      if rossiter[i] then ss.planet[i].rossiter = 1B

   endif

   ss.planet[i].fittran = fittran[i]
   ss.planet[i].fitrv = fitrv[i]
   ss.planet[i].chen = chen[i]

   ;; we can marginalize over these parameters 
   ;; even if a transit is not fit.
   ;; And with the Chen & Kipping relation, we can get a decent value
   ;; on the radius and density
   if not fittran[i] and chen[i] then ss.planet[i].cosi.scale = 1d0

   ;; Chen & Kipping prior can be used to fit planetary radius 
   ;; and we can marginalize over cosi
   if not chen[i] then begin
      if not fittran[i] then begin
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

   if i180[i] then ss.planet[i].i180 = 1


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
      printandlog, "Fitting thermal emission for " + ss.band[i].name + " band",logname
   endif

   match = where(fitreflect eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].reflect.fit = 1B
      ss.band[i].reflect.derive = 1B
      printandlog, "Fitting reflected light for " + ss.band[i].name + " band", logname
   endif

   match = where(fitdilute eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].dilute.fit = 1B
      ss.band[i].dilute.derive = 1B
      printandlog, "Fitting dilution for " + ss.band[i].name + " band", logname
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
      printandlog, 'NINTERP and EXPTIME must be unspecified, a byte, or an NTRANSITS (' + strtrim(ntran,2) + ') byte array', logname
      stop
   endelse

   printandlog, 'The index for each transit is',logname
   ss.transit[*].transitptrs = ptrarr(ntran,/allocate_heap)
   for i=0, ntran-1 do begin
      printandlog, string(i,tranfiles[i],format='(i2,x,a)'),logname
      *(ss.transit[i].transitptrs) = readtran(tranfiles[i], detrendpar=detrend)

;      stop

;      ;; create an array of detrending variables 
;      ;; (one for each extra column in the transit file)
;      nadd = (*(ss.transit[i].transitptrs)).nadd
;      if nadd ge 1 then *(ss.transit[i].detrendadd) = replicate(detrendadd,nadd)
;      nmult = (*(ss.transit[i].transitptrs)).nmult
;      if nmult ge 1 then *(ss.transit[i].detrendmult) = replicate(detrendmult,nmult)
;stop

      ss.transit[i].exptime = exptime[i]
      ss.transit[i].ninterp = ninterp[i]
      
      band = (*(ss.transit[i].transitptrs)).band
      ss.transit[i].bandndx = where(ss.band[*].name eq band)
      ss.transit[i].label = (*(ss.transit[i].transitptrs)).label
   endfor
   printandlog, '',logname
endif else begin
   ss.transit[0].f0.fit = 0
   ss.transit[0].f0.derive = 0
   ss.transit[0].variance.fit = 0
   ss.transit[0].variance.derive = 0
endelse

;; read in the RV files
if ntel gt 0 then begin
   ss.telescope[*].rvptrs = ptrarr(ntel,/allocate_heap)
   printandlog, "The index for each RV data set is",logname
   for i=0, ntel-1 do begin
      printandlog, string(i,rvfiles[i],format='(i2,x,a)'),logname
      *(ss.telescope[i].rvptrs) = readrv(rvfiles[i])
      ss.telescope[i].label = (*(ss.telescope[i].rvptrs)).label
   endfor
   printandlog, '', logname
endif else begin
   ss.telescope[*].rvptrs = ptr_new(/allocate_heap)
endelse

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
   if nentries lt 2 or nentries gt 5 then begin
      if line ne '' then printandlog, 'WARNING: line ' + strtrim(lineno,2) + ' in ' + priorfile + ' is not legal syntax (NAME VALUE [UNCERTAINTY] [LOWERBOUND] [UPPERBOUND]); ignoring: ' + line, logname
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
   if nentries ge 4 then begin
      if strupcase(entries[3]) eq '-INF' then lowerbound = -!values.d_infinity $
      else lowerbound = double(entries[3])
   endif else lowerbound = -!values.d_infinity
   if nentries ge 5 then begin
      if strupcase(entries[4]) eq 'INF' then upperbound = !values.d_infinity $
      else upperbound = double(entries[4])
   endif else upperbound = !values.d_infinity

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
               if ss.(i)[priornum].(k) ne !NULL then begin
                  for l=0L, n_tags(*(ss.(i)[priornum].(k)))-1 do begin

                     if (size((*(ss.(i)[priornum].(k))).(l)))[2] eq 8 then begin

                        valid = 0
                        if strpos(strupcase(priorlabel),'C') ne -1 then begin
                           if strmid((*(ss.(i)[priornum].(k))).(l)[0].label,0,1) eq 'C' then begin ;; if it's the additive variable
                              if (*(ss.(i)[priornum].(k))).nadd gt 0 then begin
                                 detrendnum = long((strsplit(strupcase(priorlabel),'C',/extract))[0])
                                 if detrendnum lt (*(ss.(i)[priornum].(k))).nadd then valid = 1
                              endif
                           endif
                        endif else if strpos(strupcase(priorlabel),'M') ne -1 then begin
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
                                 if nplanets ne 0 then $
                                    printandlog, "WARNING: Prior supplied on '" + $
                                                 priorname + "' but it is neither fitted or derived; not applying prior", logname
                              endif else begin
                              
                                 ;; found it! change the default starting guess to the value
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].prior = priorval
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].value = priorval
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].upperbound = upperbound
                                 (*(ss.(i)[priornum].(k))).(l)[detrendnum].lowerbound = lowerbound
                                 
                                 if priorwidth eq 0d0 then begin
                                    ;; priorwidth = 0 => fix it at the prior value
                                    (*(ss.(i)[priornum].(k))).(l)[detrendnum].fit = 0d0
                                    (*(ss.(i)[priornum].(k))).(l)[detrendnum].derive = 0d0
                                    if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (fixed)', logname
                                    
                                 endif else if finite(priorwidth) and priorwidth gt 0d0 or finite(lowerbound) or finite(upperbound) then begin
                                    ;; apply a Gaussian prior with width = priorwidth
                                    priors=[[priors],[i,priornum,k,l,detrendnum]]
                                    
                                    if priorwidth gt 0 and finite(priorwidth) then begin
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].scale = priorwidth*3d0
                                       (*(ss.(i)[priornum].(k))).(l)[detrendnum].priorwidth = priorwidth
                                    endif                  
                                    
                                    if ~keyword_set(silent) then begin
                                       if priorwidth lt 0 then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname $
                                       else printandlog, priorname + ' = ' + strtrim(priorval,2) + ' +/- ' + strtrim(priorwidth,2) + '; bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                                    endif
                                    
                                 endif else begin
                                    ;; else no prior, just change the default starting value
                                    if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
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
                  
                  if not ss.(i)[priornum].(ndx).fit and not ss.(i)[priornum].(ndx).derive then begin
                     if nplanets ne 0 then $
                        printandlog, "WARNING: Prior supplied on '" + $
                                     priorname + "' but it is neither fitted or derived; not applying prior", logname
                  endif else begin
                     
                     ;; found it! change the default starting guess to the value
                     ss.(i)[priornum].(ndx).prior = priorval
                     ss.(i)[priornum].(ndx).value = priorval
                     ss.(i)[priornum].(ndx).upperbound = upperbound
                     ss.(i)[priornum].(ndx).lowerbound = lowerbound
                     
                     if priorwidth eq 0d0 then begin
                        ;; priorwidth = 0 => fix it at the prior value
                        ss.(i)[priornum].(ndx).fit = 0d0
                        ss.(i)[priornum].(ndx).derive = 0d0
                        if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (fixed)', logname
                     endif else if finite(priorwidth) and priorwidth gt 0d0 or finite(lowerbound) or finite(upperbound) then begin
                        ;; apply a Gaussian prior with width = priorwidth
                        priors=[[priors],[i,priornum,ndx,-1,-1]]
                        
                        if priorwidth gt 0 and finite(priorwidth) then begin
                           ss.(i)[priornum].(ndx).scale = priorwidth*3d0
                           ss.(i)[priornum].(ndx).priorwidth = priorwidth
                        endif                  
                        
                        if ~keyword_set(silent) then begin
                           if priorwidth lt 0 then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname $
                           else printandlog, priorname + ' = ' + strtrim(priorval,2) + ' +/- ' + strtrim(priorwidth,2) + '; bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
                        endif
                        
                     endif else begin
                        ;; else no prior, just change the default starting value
                        if ~keyword_set(silent) then printandlog, priorname + ' = ' + strtrim(priorval,2) + ' (no prior constraint); bounded between ' + strtrim(lowerbound,2) + ' and ' +  strtrim(upperbound,2), logname
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
   if not found and nplanets ne 0 then printandlog, "WARNING: No parameter matches '" + $
      priorname + "' from " + priorfile + "; not applying prior", logname
                
endwhile

free_lun, lun
if ~keyword_set(silent) then printandlog, '', logname

;; do we have enough information to derive the distance?
;if (where(priorname eq 'distance'))[0] ne -1 
;stop

priors = priors[*,1:*]
*(ss.priors) = priors

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
            if ss.(i)[j].(k) ne !NULL then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'fit') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].fit then tofit = [[tofit],[i,j,k,l,m]]
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            ;; and this captures everything else
            if tag_exist(ss.(i)[j].(k),'fit') then begin
               if ss.(i)[j].(k).fit then tofit = [[tofit],[i,j,k,-1,-1]]
            endif
         endif
      endfor
   endfor
endfor
tofit = tofit[*,1:*]
*(ss.tofit) = tofit
ss.nchains = n_elements((*ss.tofit)[0,*])*2L

;; determine the epoch for each observation
for i=0, ntran-1 do begin
   if i eq 0 then begin
      ;; variations are defined relative to the first transit
      if keyword_set(tivs) then begin
         ss.transit[i].tiv.fit = 0
         ss.transit[i].tiv.derive = 0
      endif
      if keyword_set(tdeltavs) then begin
         ss.transit[i].tdeltav.fit = 0
         ss.transit[i].tdeltav.derive = 0
      endif
   endif

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
      
      epoch = (mean((*(ss.transit[i].transitptrs)).bjd) - tc)/period
      normepoch = ((epoch mod 1) + 1) mod 1

      ;; this causes the transit to be ignored if it's slightly
      ;; off... what was the purpose of this if statement code anyway??
;      if normepoch lt 0.05 or normepoch gt 0.95 then ss.transit[i].epoch[j] = round(epoch)

      ss.transit[i].epoch[j] = round(epoch)

   endfor
endfor

;; populate the best values
if arg_present(best) then pars2str, best, ss, /best

;; don't do these when creating the MCMC structure
if n_elements(ss.star.mstar.value) eq 1 then begin
   ;; derive all step parameters
   if not pars2step(ss) then begin
      printandlog, 'Warning: YY isochrones are not applicable here; refine priors. Are you fitting a low mass star? Be sure to disable YY isochrones using the /NOYY, disable the limb darkening prior using /NOCLARET, and supply priors on mstar, rstar, and u1 and u2 for each band', logname
      stop
   endif

   ;; return an error if the starting stellar system is not allowed
   if step2pars(ss,/verbose,/changedefaults) eq -1 then begin
      printandlog, 'Warning: starting values for stellar system not allowed; refine priors', logname
      stop
   endif

   ;; calculate the right period stepping scale for AMOEBA based on the
   ;; input data
   perscale, ss
endif

return, ss

end
