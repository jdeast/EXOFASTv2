;+
; NAME:
;   EXOFAST_CHI2
;
; PURPOSE: 
;   Computes the chi2 for a transit and/or RV for a single planet
;
; CALLING SEQUENCE:
;    chi2 = exofast_chi2(pars)
;
; INPUTS:
;
;    PARS - An array of parameters that describes the model. There
;           should be one for each parameter in the structure where
;           fit==true
; 
;      pars.rvzeros   ;; the zero points for each RV data set
;      pars.s         ;; a structure describing the star 
;        s.mass       ;; stellar mass
;        s.logg       ;; logg
;        s.teff       ;; effective temperature
;        s.feh        ;; metalicity
;        s.vsini      ;; velocity
;        s.distance   ;; stellar distance
;        s.rstar      ;; stellar radius (derived)
;        s.age        ;; stellar age (derived)
;        s.luminosity ;; stellar luminosity (derived)
;        s.umag       ;; apparent U magnitude
;        s.bmag       ;; apparent B magnitude
;        s.vmag       ;; apparent V magnitude
;        s.rmag       ;; apparent R magnitude
;        s.imag       ;; apparent I magnitude
;        s.jmag       ;; apparent J magnitude
;        s.hmag       ;; apparent H magnitude
;        s.kmag       ;; apparent K magnitude
;        s.sloanumag  ;; apparent sloan u magnitude
;        s.sloangmag  ;; apparent sloan g magnitude
;        s.sloanrmag  ;; apparent sloan r magnitude
;        s.sloanimag  ;; apparent sloan i magnitude
;        s.sloanzmag  ;; apparent sloan z magnitude
;        s.stromumag  ;; apparent stromgren u magnitude
;        s.strombmag  ;; apparent stromgren b magnitude
;        s.stromvmag  ;; apparent stromgren v magnitude
;        s.stromymag  ;; apparent stromgren y magnitude
;        s.ra         ;; right ascension
;        s.dec        ;; declination
;        s.Omega      ;; orienation in real space (??)
;        s.decpm
;        s.rapm 
;        
;     pars.p  ;; an array of structures for each planet and/or observation
;        p.tc       ;; time of inferior conjunction (primary transit)
;        p.logP     ;; log of the period
;        p.qecosw   ;; e^(1/4)*cosw
;        p.qesinw   ;; e^(1/4)*sinw
;        p.logK     ;; log of the RV semi-amplitude
;        p.cosi     ;; cosine of the inclination
;        p.p        ;; Rp/Rstar
;        p.F0       ;; Baseline flux
;        p.lambda   ;; projected Spin-Orbit alignment
;        p.dilution ;; dilution term for the planet
;        p.thermal  ;; thermal emission from the planet (constant offset)
;        p.reflect  ;; reflected flux from the planet (phase dependent)
;        p.beam     ;; beaming from the planet
;        p.ellip    ;; ellipsoidal variation of star (tidally locked to planet)
;        p.band     ;; the observed band
;        p.u1       ;; linear limb darkening (depends on band)
;        p.u2       ;; quadratic limb darkening (depends on band)
;        p.u3       ;; 1st non-linear limb darkening (depends on band)
;        p.u4       ;; 2nd non-linear limb darkening (depends on band)
;        p.slope    ;; linear trend in the RVs C*(time - t0)
;        p.quad     ;; quadratic trend in the RVs C*(time - t0)^2

;    DATA - A structure that describes the data
;    
;    data.tranptrs ;; an array pointers to each transit data set
;      (*tranptrs).trantime    ;; NDATA array of BJD_TDB
;      (*tranptrs).flux        ;; NDATA array of Normalized flux
;      (*tranptrs).fluxerr     ;; NDATA array of flux errors
;      (*tranptrs).detrendmult ;; NTRENDS x NDATA array of multiplicative trends
;      (*tranptrs).detrendsub  ;; NTRENDS x NDATA array of subtractive trends
;    data.rvptrs  ;; an array of pointers to each RV data set
;    data.options ;; a structure of options on how to fit the data
;      options.fitrvs  ;; an NRV array of which RV data sets to fit  
;      options.fittran ;; an NPLANET array of which transit data sets to fit  
;      options.band    ;; an NPLANET array of the observed bands
;      options.exptime ;; an NPLANET array of exposure times
;      options.ninterp ;; an NPLANET array of number of interpolations
;      options.circular ;; an NPLANET array of which planets should be circular
;      options.yy      ;; use the YY isochrones
;      options.tides   ;; include a rough treatment of tidal circularization
;      options.debug   ;; debug

;    PRIORS - A structure that mirrors PARS but with two-element arrays
;             for each parameter [prior, prior width].
;    LABELS - A structure that mirrors PARS but with 3-element strings
;             arrays for each parameter. Each parameter is a 3-element array:
;             Human-readable label
;             Human-readable explanation
;             textoidl form for machine-readable plotting
;             labels.gamma[0] = 'gamma'
;             labels.gamma[1] = 'Systemic velocity (m/s)'
;             labels.gamma[2] = textoidl('\gamma')

;    PARS - a parameter array containing all of the parameters in the
;           model.
;
;           gamma     = pars[0]       ;; systemic velocity
;           slope     = pars[1]       ;; slope in RVs
;           tc        = pars[2]       ;; transit center time
;           logP      = pars[3]       ;; alog10(Period/days)
;           qecosw    = pars[4]       ;; eccentricity/arg of periastron
;           qesinw    = pars[5]       ;; eccentricity/arg of periastron
;           logK      = pars[6]       ;; alog10(velocity semi-amplitude/(m/s))
;           cosi      = pars[7]       ;; cosine of inclination of the orbit
;           p         = pars[8]       ;; rp/rstar
;           log(ar)   = pars[9]       ;; alog10(a/rstar)
;           logg      = pars[10]      ;; stellar surface gravity
;           teff      = pars[11]      ;; stellar effective temperature
;           feh       = pars[12]      ;; stellar metallicity
;           depth2    = pars[13]      ;; secondary eclipse depth
;           u1        = pars[14]      ;; linear limb darkening coeff
;           u2        = pars[15]      ;; quadratic limb darkening coeff
;           u3        = pars[16]      ;; 1st non-linear limb darkening coeff (not supported)
;           u4        = pars[17]      ;; 2nd non-linear limb darkening coeff (not supported)
;           F0        = pars[18]      ;; baseline flux
;           coeffs = pars[19:npars-1] ;; detrending variables
;
; OPTIONAL INPUTS:
;    PSNAME      - The name of a PS file. If set, a plot the
;                  data/model will be written to this file.
; OPTIONAL OUTPUTS:
;    DETERMINANT - The determinant of the parameterization above and
;                  the uniform priors we wish to impose. In this case,
;                  it is always 1d0 (but is required by EXOFAST_DEMC).
;    MODELRV     - The RV model at each time (rv.bjd).
;    MODELFLUX   - The model light curve at each time (transit.bjd).
;   
;
; RESULT:
;    The chi^2 of the model given the data and parameters.
;
; COMMON BLOCKS:
;   CHI2_BLOCK - See exofast.pro for definition
;
; MODIFICATION HISTORY
; 
;  2012/06 -- Public release -- Jason Eastman (LCOGT)
;  2012/07 -- Fixed major bug in mstar/rstar prior width derivation
;  2012/12 -- Add Long cadence, quadratic limb darkening fit.
;  2012/12 -- Changed eccentricity constraint to e < (1-Rstar/a)
;  2013/02 -- Fixed bug that broke detrending, introduced in 2012/12
;  2013/03 -- Changed eccentricity prior:
;               e > (1-Rstar/a-rp/a) -- models are rejected
;               e > (1-3*Rstar/a) -- eccentricity set to zero if options.tides=1
;               now step in e^(1/4)*cos(omega), e^(1/4)*sin(omega) to
;               more closely match the observed prior distribution.
;             Added support for Mstar/Rstar priors, replaces Torres relation
;             Added support for YY evolutionary models, replaces
;             Torres relation
;             Added DERIVED parameter, which returns the age.
;-

function exofast_chi2v2, pars, determinant=determinant, $
                         modelrv=modelrv, modelflux=modelflux, psname=psname, $
                         derived=derived
  
COMMON chi2_block, ss
;; populate the structure with the new parameters
if n_elements(pars) ne 0 then pars2str, pars, ss

au = ss.constants.au/ss.constants.rsun

;; derive all required parameters 
;; (this may change depending on parameterization)
;derivepars, ss

;; initialize the determinant and chi^2
chi2 = 0.d0
determinant = 1d0

;; physical limb darkening
if (where(ss.planet.fittran))[0] ne -1 then begin
   bad = where(ss.band.u1.value + ss.band.u2.value gt 1d0 or $
               ss.band.u1.value + ss.band.u2.value lt 0d0 or $
               ss.band.u2.value lt -1 or ss.band.u2.value gt 1d0, nbad)
   
   if nbad gt 0 then begin
      if ss.debug then printandlog, strtrim(nbad,2) + ' limb darkening parameters are bad (' + strtrim(ss.band[bad[0]].u1.value,2) + ', ' + strtrim(ss.band[bad[0]].u2.value,2) + ')',ss.logname
      return, !values.d_infinity
   endif
endif

;; prevent runaways
bad = where(ss.planet.logp.value gt 7d0 or ss.planet.logp.value lt -1d0,nbad)
;bad = where(ss.planet.logp.value gt 3d0 or ss.planet.logp.value lt 0d0,nbad)
if nbad gt 0 then begin
   if ss.debug then printandlog, 'logP is bad (' + strtrim(bad,2) + ')', ss.logname
   return, !values.d_infinity
endif

;; prevent runaways
bad = where(ss.planet.tc.value lt ss.planet.tc.prior - 10^ss.planet.logp.value or $
            ss.planet.tc.value gt ss.planet.tc.prior + 10^ss.planet.logp.value,nbad)
if nbad gt 0 then begin
   if ss.debug then printandlog, 'tc is bad (' + strtrim(bad,2) + ')', ss.logname
   return, !values.d_infinity
endif

;; 0 <= cosi <= 1 (or -1 <= cosi <= 1 if i180 keyword set)
bad = where(ss.planet.cosi.value gt 1 or (ss.planet.cosi.value lt 0 and ~ss.planet.i180) or (ss.planet.cosi.value lt -1),nbad)
if nbad gt 0 then begin
   if ss.debug then printandlog, 'cosi is bad (' + strtrim(bad,2) + ')', ss.logname
   return, !values.d_infinity
endif

;; older than the universe (too conservative?)
if ss.star.age.value gt 13.82d0 or ss.star.age.value lt 0d0 then begin
   if ss.debug then printandlog, strtrim(nbad,2) + ' age is bad (' + strtrim(bad,2) + ')',ss.logname
   return, !values.d_infinity
endif

;; positive extinction
if ss.star.av.value lt 0 then begin
   if ss.debug then printandlog, 'extinction is bad (' + strtrim(ss.star.av.value,2) + ')',ss.logname
   return, !values.d_infinity
endif

;; positive distance
if ss.star.distance.value lt 0 then begin
   if ss.debug then printandlog, 'distance is bad (' + strtrim(ss.star.distance.value,2) + ')',ss.logname
   return, !values.d_infinity
endif

;; bound marginally detected planets to limit (infinite) parameter space at low logK
;; conservative lower limit of 0.1 Ceres in 1 year orbit around sun = 1 um/s
;; conservative upper limit corresponds to ~thousand solar masses
bad = where((ss.planet.logk.value lt -6d0 or ss.planet.logk.value gt 6d0) and ss.planet.logk.fit, nbad)
if nbad gt 0 then begin
   if ss.debug then printandlog, 'k is bad (' + strtrim(bad,2) + ')', ss.logname
   return, !values.d_infinity
endif

;; derive the model parameters from the stepping parameters (return if unphysical)
ss.star.mstar.value = 10^ss.star.logmstar.value

;; use the YY tracks to guide the stellar parameters
if ~ss.noyy then begin
   if keyword_set(psname) then epsname = psname+'.yy.eps'
   yychi2 = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, $
                           ss.star.age.value, ss.star.teff.value,$
                           yyrstar=yyrstar, debug=ss.debug, psname=epsname, $
                           sigmab=ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2, $
                           gravitysun=ss.constants.gravitysun) 
   if ~finite(yychi2) then begin
      if ss.debug then printandlog, 'star is bad', ss.logname
      return, !values.d_infinity
   endif

   yychi2 += ((ss.star.rstar.value - yyrstar)/(0.03*yyrstar))^2
   chi2 += yychi2
   if ss.verbose then printandlog, 'YY penalty = ' + strtrim(yychi2,2), ss.logname
endif

;; use the MIST tracks to guide the stellar parameters
if ss.mist then begin
   if keyword_set(psname) then epsname = psname+'.mist.eps'
   mistchi2 = massradius_mist(ss.star.eep.value, ss.star.mstar.value, ss.star.initfeh.value, $
                              ss.star.age.value, ss.star.teff.value,$
                              ss.star.rstar.value, ss.star.feh.value, debug=ss.debug, $
                              epsname=epsname, gravitysun=ss.constants.gravitysun)
   chi2 += mistchi2
   if ss.verbose then printandlog, 'MIST penalty = ' + strtrim(mistchi2,2), ss.logname
   if ~finite(chi2) then begin
      if ss.debug then printandlog, 'star not on MIST tracks', ss.logname
      return, !values.d_infinity
   endif
endif

;; use the Torres relation to guide the stellar parameters
if ss.torres then begin
   massradius_torres, ss.star.logg.value, ss.star.teff.value, ss.star.feh.value, mstar_prior, rstar_prior
   umstar = 0.027d0;*100d0 ;; why the factor of 100?!
   urstar = 0.014d0;*100d0
   if mstar_prior lt 0.6d0 then printandlog, $
      'WARNING: Torres not applicable (mstar = ' + $
      strtrim(mstar_prior,2) + '); ignore at beginning. Otherwise, ' + $
      'use MIST, YY, or impose a prior on mstar/rstar',ss.logname
   ;; add "prior" penalty
   chi2 += (alog10(ss.star.mstar.value/mstar_prior)/umstar)^2
   chi2 += (alog10(ss.star.rstar.value/rstar_prior)/urstar)^2

   if ss.verbose then $
      printandlog, 'Torres penalty: ' + string((alog10(ss.star.mstar.value/mstar_prior)/umstar)^2,$
                                               (alog10(ss.star.rstar.value/rstar_prior)/urstar)^2), ss.logname

endif

if ss.star.errscale.value le 0 then begin
   if ss.debug then printandlog, 'error scale is bad', ss.logname
   return, !values.d_infinity
endif

bad = where(ss.doptom.dtscale.value le 0, nbad)
if nbad gt 0 then begin
   if ss.debug then printandlog, 'dtscale is bad (' + strtrim(bad,2) + ')', ss.logname
   return, !values.d_infinity
endif

if ss.star.alpha.value lt -0.3d0 or ss.star.alpha.value gt 0.7d0 then begin
   if ss.debug then printandlog, 'alpha is bad (' + strtrim(ss.star.alpha.value,2) + ')', ss.logname
   return, !values.d_infinity
endif

if step2pars(ss,verbose=ss.debug,logname=ss.logname) eq -1 then begin
   if ss.debug then printandlog, 'stellar system is bad', ss.logname
   return, !values.d_infinity
endif

;; add prior penalties
priors = *(ss.priors)
for i=0, n_elements(priors[0,*])-1 do begin

   ;; apply user-defined bounds
   if ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value gt ss.(priors[0,i])[priors[1,i]].(priors[2,i]).upperbound or $
      ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value lt ss.(priors[0,i])[priors[1,i]].(priors[2,i]).lowerbound then begin
      if ss.debug then $
         printandlog, ss.(priors[0,i])[priors[1,i]].(priors[2,i]).label + '_' + strtrim(priors[1,i],2) + $
                      ' (' + strtrim(ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value,2) + ') is out of bounds (' +$
                      strtrim(ss.(priors[0,i])[priors[1,i]].(priors[2,i]).upperbound,2) + ',' + strtrim(ss.(priors[0,i])[priors[1,i]].(priors[2,i]).lowerbound,2) + ')',ss.logname
      return, !values.d_infinity
   endif

   
   chi2 += ((ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value - $
             ss.(priors[0,i])[priors[1,i]].(priors[2,i]).prior)/$
            ss.(priors[0,i])[priors[1,i]].(priors[2,i]).priorwidth)^2
   
;   printandlog, ss.(priors[0,i])[priors[1,i]].(priors[2,i]).label,  ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value, $
;          ss.(priors[0,i])[priors[1,i]].(priors[2,i]).prior, $
;          ss.(priors[0,i])[priors[1,i]].(priors[2,i]).priorwidth, $
;          ((ss.(priors[0,i])[priors[1,i]].(priors[2,i]).value - $
;            ss.(priors[0,i])[priors[1,i]].(priors[2,i]).prior)/$
;           ss.(priors[0,i])[priors[1,i]].(priors[2,i]).priorwidth)^2,$
;          ss.(priors[0,i])[priors[1,i]].(priors[2,i]).lowerbound,$ 
;          ss.(priors[0,i])[priors[1,i]].(priors[2,i]).upperbound, ss.logname

endfor

;; Apply the Mass-Radius relation 
;; Chen & Kipping, 2017 (http://adsabs.harvard.edu/abs/2017ApJ...834...17C)
;; this introduces a near perfect correlation between K and p, which
;; AMOEBA finds challenging to work with (but is handled naturally by DEMC)
;; does it? Or was that a bug that has been fixed?
for j=0, ss.nplanets-1 do begin
   if ss.planet[j].chen then begin

      ;; negative radii are allowed to assess the significance of the
      ;; transit depth. That breaks these relations, so exclude them here
      if ss.planet[j].rpearth.value le 0d0 then begin
         if ss.debug then printandlog, 'rpearth is bad', ss.logname
         return, !values.d_infinity
      endif

      rp = massradius_chen(ss.planet[j].mpearth.value,rperr=rperr)

      ;; add a chi2 penalty for deviation from the mass-radius relation
      ;; if the radius is well-constrained (by transit depth), it
      ;; becomes an implicit constraint on mass. If the mass is well
      ;; constrained (by RV), it becomes an implicit constraint on
      ;; radius
      chi2 += ((rp - ss.planet[j].rpearth.value)/rperr)^2
      
;      printandlog, 'chen penalty = ' + strtrim(((rp - ss.planet[j].rpearth.value)/rperr)^2,2),logname
;      printandlog, rp, ss.planet[j].rpearth.value,rperr,logname
;      printandlog, ss.planet[j].mpearth.value,logname
;      printandlog, ss.planet[j].k.value,logname

   endif
endfor

;; fit the SED
if file_test(ss.star.fluxfile) then begin
   if keyword_set(psname) then epsname = psname+'.sed.eps'
   sedchi2 = exofast_sed(ss.star.fluxfile, ss.star.teff.value, $
                         ss.star.rstar.value,$
                         ss.star.av.value, ss.star.distance.value, $
                         logg=ss.star.logg.value,met=ss.star.feh.value,$
                         alpha=ss.star.alpha.value,verbose=ss.debug, $
                         f0=f, fp0=fp, ep0=ep, psname=epsname, $
                         pc=ss.constants.pc, rsun=ss.constants.rsun)

   if ~finite(sedchi2) then begin
      if ss.debug then printandlog, 'sed is bad', ss.logname
      return, !values.d_infinity
   endif

   sedchi2 = exofast_like(f-fp,0d0,ss.star.errscale.value*ep,/chi2)
   if ~finite(sedchi2) then begin
      if ss.debug then printandlog, 'sed is bad', ss.logname
      return, !values.d_infinity
   endif
   chi2 += sedchi2
;   printandlog, 'SED penalty = ' + strtrim(sedchi2,2), ss.logname
endif

;; RV model (non-interacting planets)
for j=0, ss.ntel-1 do begin

   rv = *(ss.telescope[j].rvptrs)

   if (where(rv.err^2 + ss.telescope[j].jittervar.value le 0d0))[0] ne -1 then return, !values.d_infinity

   modelrv = dblarr(n_elements(rv.rv))
   for i=0, ss.nplanets-1 do begin

      if ss.planet[i].fitrv then begin      
         ;; rvbjd = rv.bjd ;; usually sufficient (See Eastman et al., 2013)

         ;; time in target barycentric frame (expensive)
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value, $
                            a=ss.planet[i].a.value, tp=ss.planet[i].tp.value, $
                            period=ss.planet[i].period.value, e=ss.planet[i].e.value,$
                            omega=ss.planet[i].omega.value,/primary)
         
         ;; calculate the RV model
         modelrv += exofast_rv(rvbjd,ss.planet[i].tp.value,ss.planet[i].period.value,$
                               0d0,ss.planet[i].K.value,$
                               ss.planet[i].e.value,ss.planet[i].omega.value,$
                               slope=0, $
                               rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value,a=ss.planet[i].ar.value,$
                               p=ss.planet[i].p.value,vsini=ss.star.vsini.value,$
                               lambda=ss.planet[i].lambda.value,$
                               u1=0d0,t0=t0,deltarv=deltarv)

      endif

   endfor
   ;; add instrumental offset, slope, and quadratic term
   modelrv += ss.telescope[j].gamma.value + ss.star.slope.value*(rv.bjd-t0) + ss.star.quad.value*(rv.bjd-t0)^2

   (*ss.telescope[j].rvptrs).residuals = rv.rv - modelrv
   
   if keyword_set(psname) then $
      exofast_forprint, rv.bjd, rv.rv - modelrv, rv.err, format='(f0.8,x,f0.6,x,f0.6)', textout=psname + '.residuals.' + strtrim(j,2) + '.rv.txt', /nocomment,/silent
   if keyword_set(psname) then $
      exofast_forprint, rv.bjd, modelrv, format='(f0.8,x,f0.6)', textout=psname + '.model.' + strtrim(j,2) + '.rv.txt', /nocomment,/silent

;   printandlog, ss.telescope[j].jittervar.value,ss.telescope[j].jitter.value, sqrt(ss.telescope[j].jittervar.value),ss.logname
   rvchi2 = exofast_like((*ss.telescope[j].rvptrs).residuals,ss.telescope[j].jittervar.value,rv.err,/chi2)
   
   if ~finite(rvchi2) then stop
   chi2 += rvchi2
;   printandlog, 'RV penalty = ' + strtrim(rvchi2,2),ss.logname
;   chi2 += total(((rv.rv - modelrv)/rv.err)^2)
endfor

;; if at least one RV planet is fit, plot it
if (where(ss.planet.fitrv))[0] ne -1 then begin
   if keyword_set(psname) then begin
      plotrv, ss, psname=psname + '.rv.ps'
   endif else if ss.debug then begin
      plotrv, ss
   endif
endif

;; Doppler Tomography Model
for i=0, ss.ndt-1 do begin
   if keyword_set(psname) then epsname = psname + 'dt.eps'
   dtchi2 = dopptom_chi2(*(ss.doptom[i].dtptrs),$
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].tc.value, $
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].period.value, $
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].e.value,$
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].omega.value, $
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].cosi.value, $
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].p.value,$
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].ar.value,$
                         ss.planet[(*(ss.doptom[i].dtptrs)).planetndx].lambda.value, $
                         ss.star.logg.value, ss.star.teff.value, ss.star.feh.value,$
                         ss.star.vsini.value/1d3,ss.star.macturb.value/1d3,$
                         ss.doptom[i].dtscale.value, debug=ss.debug,/like,psname=epsname)

   if ~finite(dtchi2) then return, !values.d_infinity
   chi2 += dtchi2
endfor

;; Transit model
for j=0, ss.ntran-1 do begin

   transit = *(ss.transit[j].transitptrs)

   if (where(transit.err^2 + ss.transit[j].variance.value le 0d0))[0] ne -1 then return, !values.d_infinity

   band = ss.band[ss.transit[j].bandndx]

   ;; quadratic limb darkening
   if ~ss.noclaret then begin
      ldcoeffs = quadld(ss.star.logg.value, ss.star.teff.value, ss.star.feh.value, band.name)
      u1claret = ldcoeffs[0]
      u2claret = ldcoeffs[1]
      u1err = 0.05d0 
      u2err = 0.05d0
      chi2 += ((band.u1.value-u1claret)/u1err)^2
      chi2 += ((band.u2.value-u2claret)/u2err)^2
;      printandlog, 'u1 penalty = ' + strtrim(((band.u1.value-u1claret)/u1err)^2,2),ss.logname
;      printandlog, 'u2 penalty = ' + strtrim(((band.u2.value-u2claret)/u2err)^2,2),ss.logname
   endif

   ;; Kepler Long candence data; create several model points and average   
   ninterp = ss.transit[j].ninterp
   npoints = n_elements(transit.bjd)
   if ninterp gt 1 then begin
      transitbjd = transit.bjd#(dblarr(ninterp)+1d0) + $     
                   ((dindgen(ninterp)/(ninterp-1d0)-0.5d0)/$
                    1440d*ss.transit[j].exptime)##(dblarr(npoints)+1d)
      modelflux = dblarr(npoints,ninterp) + 1d0
   endif else begin
      transitbjd = transit.bjd
      modelflux = dblarr(npoints) + 1d0
   endelse

   ;; get the motion of the star due to the planet
   junk = exofast_getb2(transitbjd,inc=ss.planet.i.value,a=ss.planet.ar.value,$
                        tperiastron=ss.planet.tp.value,$
                        period=ss.planet.period.value,$
                        e=ss.planet.e.value,omega=ss.planet.omega.value,$
                        q=ss.star.mstar.value/ss.planet.mpsun.value,$
                        x1=x1,y1=y1,z1=z1)

   for i=0, ss.nplanets-1 do begin
      if ss.planet[i].fittran then begin
         
;         printandlog, ss.planet[0].p.value, ss.logname

         modelflux += (exofast_tran(transitbjd, $
                                    ss.planet[i].i.value + ss.transit[j].tiv.value, $
                                    ss.planet[i].ar.value, $
                                    ss.planet[i].tp.value + ss.transit[j].ttv.value, $
                                    ss.planet[i].period.value, $
                                    ss.planet[i].e.value,$
                                    ss.planet[i].omega.value,$
                                    ss.planet[i].p.value + ss.transit[j].tdeltav.value,$
                                    band.u1.value, $
                                    band.u2.value, $
                                    1d0, $
                                    q=ss.star.mstar.value/ss.planet[i].mpsun.value, $
                                    thermal=band.thermal.value, $
                                    reflect=band.reflect.value, $
                                    dilute=band.dilute.value,$
                                    tc=ss.planet[i].tc.value,$
                                    rstar=ss.star.rstar.value/AU,$
                                    x1=x1,y1=y1,z1=z1,au=au) - 1d0)
      endif

   endfor

   modelflux *=  ss.transit[j].f0.value
   
   ;; now integrate the model points (before detrending)
   if ninterp gt 1 then modelflux = total(modelflux,2)/ninterp

;   if band.thermal.fit then begin
;      plot, transitbjd, transit.flux, /ynoz
;      oplot, transitbjd, modelflux,color='0000ff'x
;   endif

   ;; detrending
   modelflux += total(transit.detrendadd*replicate(1d0,n_elements(transit.bjd))##transit.detrendpars.value,1)  

   ;; chi^2
   transitchi2 = exofast_like(transit.flux - modelflux,ss.transit[j].variance.value,transit.err,/chi2)


   (*ss.transit[j].transitptrs).residuals = transit.flux - modelflux
   (*ss.transit[j].transitptrs).model = modelflux

   if keyword_set(psname) then $
      exofast_forprint, transit.bjd, transit.flux - modelflux, transit.err, format='(f0.8,x,f0.6,x,f0.6)', textout=psname + '.residuals.' + strtrim(j,2) + '.flux.txt', /nocomment,/silent
   if keyword_set(psname) then $
      exofast_forprint, transit.bjd, modelflux, format='(f0.8,x,f0.6)', textout=psname + '.model.' + strtrim(j,2) + '.flux.txt', /nocomment,/silent

   if ~finite(transitchi2) then stop

;transitchi2 = strtrim(total(((transit.flux - modelflux)/transit.err)^2),2)

   chi2 += transitchi2
;   printandlog, 'transit penalty: ' + strtrim(transitchi2,2) + ' ' + strtrim(ss.transit[j].variance.value,2),ss.logname

;printandlog, 'transit chi2 = ' + strtrim(total(((transit.flux - modelflux)/transit.err)^2),2),ss.logname
;stop
;   screen = GET_SCREEN_SIZE()
;   if win_state(5) then wset, 5 $
;   else window, 5, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
;   transitchi22 = total((transit.flux - modelflux)^2/(ss.transit[j].variance.value+transit.err^2))
;   plot, transitbjd, transit.flux-modelflux, title=strtrim(transitchi2,2) +  ' ' + strtrim(transitchi22,2), yrange=[-0.001,0.001]


endfor

;; plot the transit model and data 
;; if a transit is fit for at least one planet
if ((where(ss.planet.fittran))[0] ne -1) then begin
   if keyword_set(psname) then begin
      plottran, ss, psname=psname + '.transit.ps'
   endif else if ss.debug then begin
      plottran, ss
   endif
endif

;; if TTVs are allowed, add a chi2 penalty to the period and t0 from
;; the fit to a linear ephemeris
;; this imposes the constraint of the linear ephemeris while allowing TTVs
if ss.ttvs then begin
   ;; an nplanets x ntransits array of model transit times
   time = replicate(1,ss.nplanets)#ss.transit.ttv.value + $
          ss.planet.tc.value#replicate(1,ss.ntran) + $
          ss.planet.period.value#replicate(1,ss.ntran)*$
          ss.transit.epoch
   ;; add a chi2 penalty for the deviation of the ephemeris to the
   ;; linear fit of the transit times for each planet
   for i=0L, ss.nplanets-1L do begin
      good = where(finite(time[i,*]),ngood)
      if ngood lt 2 then continue

      coeffs = poly_fit((ss.transit.epoch)[i,good],time[i,good],1, sigma=sigma, yfit=yfit)
      sigma = (sigma > 1d-18)

      chi2 += ((coeffs[0]-ss.planet[i].tc.value)/sigma[0])^2
      chi2 += ((coeffs[1]-ss.planet[i].period.value)/sigma[1])^2
      
;      printandlog, ((coeffs[0]-ss.planet[i].tc.value)/sigma[0])^2,ss.logname
;      printandlog, ((coeffs[0]-ss.planet[i].tc.value)/sigma[0])^2,ss.logname

      if ss.debug or keyword_set(psname) then begin
         if keyword_set(psname) then begin
            ;; astrobetter.com tip on making pretty IDL plots
            mydevice=!d.name
            set_plot, 'PS'
            aspect_ratio=1.5
            xsize=10.5
            ysize=xsize/aspect_ratio
            !p.font=0
            device, filename=psname + '.ttv.' + strtrim(i,2) + '.eps', /color, bits=24,/encapsulated
            device, xsize=xsize, ysize=ysize
            LOADCT, 39,/silent
            colors = [0,254,159,95,223,31,207,111,191,47]
            charsizelegend = 0.09
            xlegend = 0.1
            ylegend = 0.90
            charsize = 0.5
         endif else begin
            device,window_state=win_state
            if win_state[10+i] eq 1 then wset, 10+i $
            else window, 10+i, retain=2
            colors= ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x,'0080ff'x,$
                     '800080'x,'00ffff'x,'ffff00'x,'80d000'x,'660000'x]
            charsizelegend = 0.03
            xlegend = 0.90
            ylegend = 0.95
            charsize = 1
         endelse
         ncolors = n_elements(colors)
         syms = [0,3,8,5,0,3,8,5]
         fill = [1,1,1,1,0,0,0,0]
         nsyms = n_elements(syms)
         
         telescopes = strarr(ss.ntran)
         for j=0L, ss.ntran-1L do telescopes[j] = (strsplit(ss.transit[j].label,' UT ',/regex,/extract))[0]
         sorted = sort(telescopes)
         tnames = telescopes[sorted[uniq(telescopes[sorted])]]

         xmin = min((ss.transit.epoch)[i,good],max=xmax)
         ymin = min((time[i,good]-yfit)*86400d0,max=ymax)
         plot, [0],[0],psym=3, xtitle='!3Epoch', ytitle='!3O-C (seconds)',xrange=[xmin,xmax],yrange=[ymin,ymax]
         for j=0, n_elements(tnames)-1 do begin
            observed = where(telescopes eq tnames[j])
            if observed[0] ne -1 then begin
               plotsym, syms[j mod nsyms], color=colors[j mod ncolors],fill=fill[j mod nsyms]
               oplot, ((ss.transit.epoch)[i,good])[observed],(time[i,good[observed]]-yfit[observed])*86400d0,psym=8
               xsize = (!x.crange[1] - !x.crange[0])
               ysize = (!y.crange[1] - !y.crange[0])

               ;; only need a legend if we have more than one telescope
               if n_elements(tnames) gt 1 then begin
                  xyouts, !x.crange[0] + xlegend*xsize,!y.crange[0]+(ylegend - j*charsizelegend)*ysize, $
                          tnames[j],color=colors[j mod ncolors],charsize=charsize
                  oplot, [!x.crange[0]+xlegend*xsize-xsize/20],$
                         [!y.crange[0]+(ylegend - (j-0.25)*charsizelegend)*ysize],psym=8
               endif
            endif
         endfor
         oplot, [-9d9,9d9],[0d0,0d0],linestyle=2

         if keyword_set(psname) then begin
            !p.font=-1
            !p.multi=0
            device, /close
            device, encapsulated=0
            set_plot, mydevice
         endif

;         printandlog, ((coeffs[0]-ss.planet[i].tc.value)/sigma[0])^2, ((coeffs[1]-ss.planet[i].period.value)/sigma[1])^2,ss.logname
      endif
   endfor
endif


;   plot, transitbjd, transit.flux, psym=1,/ynoz
;   oplot, transitbjd, modelflux, color=red

;; print all the parameters and the chi^2
if ss.debug then printandlog, string(pars, chi2, format='(' + strtrim(n_elements(pars)+1,2) + '(f0.8,x))'),ss.logname

;printandlog, ss.planet[0].p.value, ss.logname
;wait, 0.1
;stop

;if keyword_set(psname) then begin
;   device, /close
;   set_plot, mydevice
;endif

;; if this stop is triggered, you've found a bug!!
if ~finite(chi2) then stop

return, chi2

end

