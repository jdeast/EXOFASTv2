;; this derives all relevant parameters from the fitted parameters in
;; the stellar structure
function step2pars, ss, verbose=verbose, logname=logname, changedefaults=changedefaults

G = ss.constants.GMSun/ss.constants.RSun^3*ss.constants.day^2 ;; R_sun^3/(m_sun*day^2)
AU = ss.constants.au/ss.constants.rsun ;; R_sun
mjup = ss.constants.gmjupiter/ss.constants.gmsun ;; m_sun
rjup = ss.constants.rjupiter/ss.constants.rsun  ;; r_sun
mearth = ss.constants.gmearth/ss.constants.gmsun ;; m_sun
rearth = ss.constants.rearth/ss.constants.rsun  ;; r_sun
sigmaB = ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; Stefan-Boltzmann constant

;;; these are the default guesses... why are these here?
;if keyword_set(changedefaults) then begin
;   if ss.star.rstar.value ne 1d0 and ss.star.rstarsed.value eq 1d0 then ss.star.rstarsed.value = ss.star.rstar.value 
;   if ss.star.rstar.value eq 1d0 and ss.star.rstarsed.value ne 1d0 then ss.star.rstar.value = ss.star.rstarsed.value 
;   
;   if ss.star.teff.value ne 5778d0 and ss.star.teffsed.value eq 5778d0 then ss.star.teffsed.value = ss.star.teff.value 
;   if ss.star.teff.value eq 5778d0 and ss.star.teffsed.value ne 5778d0 then ss.star.teff.value = ss.star.teffsed.value 
;endif

;; derive stellar parameters
ss.star.mstar.value = 10^ss.star.logmstar.value
ss.star.logg.value = alog10(ss.constants.gravitysun*ss.star.mstar.value/(ss.star.rstar.value^2)) ;; cgs
;; derive the distance from lstar
ss.star.lstar.value = 4d0*!dpi*ss.star.rstar.value^2*ss.star.teff.value^4*sigmaB ;; L_sun
if ss.star.distance.fit then ss.star.parallax.value = 1d3/ss.star.distance.value ;; mas
if ss.star.parallax.fit then ss.star.distance.value = 1d3/ss.star.parallax.value ;; mas
ss.star.fbol.value = (ss.star.lstar.value/ss.constants.lsun)/(4d0*!dpi*(ss.star.distance.value/ss.constants.pc)^2) ;; cgs
ss.star.rhostar.value = ss.star.mstar.value/(ss.star.rstar.value^3)*ss.constants.rhosun ;; rho_sun

for j=0, ss.ntel-1 do begin
   if ss.telescope[j].jittervar.value gt 0 then $
      ss.telescope[j].jitter.value = sqrt(ss.telescope[j].jittervar.value)
endfor

for i=0, ss.nplanets-1 do begin

   ;; derive the mass of the planet
   if ss.planet[i].logmp.fit then ss.planet[i].mpsun.value = 10^ss.planet[i].logmp.value $ ;; m_sun
   else ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value) ;; m_sun
   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup ;; m_jupiter
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth ;; m_earth

   ;; derive the radius of the planet   
   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star.rstar.value ;; r_sun
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup ;; r_jupiter
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth ;; r_earth

   ;; allow neg planet mass to avoid boundary bias, but can't have neg total mass
   if ss.star.mstar.value + ss.planet[i].mpsun.value le 0d0 then begin
      if keyword_set(verbose) then printandlog, 'Net mass (Mstar + Mp) must be positive; rejecting step', logname
      return, -1
   endif

   ;; derive period
   ss.planet[i].period.value = 10^ss.planet[i].logp.value

   ;; possible e/omega/i/mp parameterizations are:
   ;; secosw, sesinw, cosi
   ;; ecosw, esinw
   ;; vcve, lsinw, lcosw, cosi/chord, sign
   ;; vcve, lsinw, lcosw, cosi/chord
   ;; e, omega, cosi

   ;; derive eccentricity and argument of periastron
   if ss.planet[i].secosw.fit and ss.planet[i].sesinw.fit then begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endif else if ss.planet[i].vcve.fit then begin

      if ss.planet[i].lsinw2.fit then begin
         if (ss.planet[i].vcve.value-1d0 lt 0) then $
            ss.planet.lsinw.value = -ss.planet.lsinw2.value $
         else ss.planet.lsinw.value = ss.planet.lsinw2.value
      endif

      ;; based on L*cos(omega) and L*sin(omega)
      if (ss.planet[i].lsinw2.fit or ss.planet[i].lsinw.fit) and ss.planet[i].lcosw.fit then begin
         ss.planet[i].omega.value = atan(ss.planet[i].lsinw.value, ss.planet[i].lcosw.value)
      endif

      ;; what sign of the quadratic solution?
      if ss.planet[i].sign.fit then begin
         ;; we fit for it
         sign = floor(ss.planet[i].sign.value)
         ;; ****** this flips the sign of the sign and is experimental*******
         if ss.planet[i].vcve.value lt 1d0 then sign = ~sign
      endif else begin
         ;; L implicitly defines it
         lsq = ss.planet[i].lsinw.value^2 + ss.planet[i].lcosw.value^2
         sign = (lsq ge 0.5d0) ;; if L^2 >= 0.5, use the positive solution to the quadratic!
      endelse

      ss.planet[i].e.value = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=sign)
      if ~finite(ss.planet[i].e.value) or ss.planet[i].e.value lt 0d0 or ss.planet[i].e.value ge 1d0 then begin
         if keyword_set(verbose) then printandlog, 'e is not in range', logname
         return, -1   
      endif
   endif else if ss.planet[i].ecosw.fit and ss.planet[i].esinw.fit then begin
      ;; ecosw, esinw parameterization
      ss.planet[i].e.value = sqrt(ss.planet[i].ecosw.value^2 + ss.planet[i].esinw.value^2)
      ss.planet[i].omega.value = atan(ss.planet[i].esinw.value,ss.planet[i].ecosw.value)
   endif else if ss.planet[i].e.fit and ss.planet[i].omega.fit then begin
      ;; fit e, omega direct
      if ss.planet[i].e.value lt 0d0 or ss.planet[i].e.value ge 1d0 then begin
         if keyword_set(verbose) then printandlog, 'Eccentricity not allowed'
         return, -1
      endif
      ;; scale omega from -pi to pi
      ss.planet[i].omega.value = (ss.planet[i].omega.value mod (2d0*!dpi)) 
      if ss.planet[i].omega.value gt !dpi then ss.planet[i].omega.value -= 2d0*!dpi
      if ss.planet[i].omega.value le -!dpi then ss.planet[i].omega.value += 2d0*!dpi
   endif else if ss.planet[i].qecosw.fit and ss.planet[i].qesinw.fit then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
      ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)
   endif

   ;; derive various combinations of e and omega
   if ss.planet[i].e.value eq 0d0 then ss.planet[i].omega.value = !dpi/2d0
   ss.planet[i].omegadeg.value = ss.planet[i].omega.value*180d0/!dpi
   if ~ss.planet[i].esinw.fit then ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].ecosw.fit then ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].vcve.fit then ss.planet[i].vcve.value = sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; derive a/rstar, a
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                       ;; semi-major axis in r_sun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value    ;; a/rstar (unitless)
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; semi major axis in AU

   ;; limit eccentricity to avoid collision with star during periastron
   ;; the ignored tidal effects would become important long before this,
   ;; but this prevents numerical problems compared to the e < 1 constraint
   ;; the not/lt (instead of ge) robustly handles NaNs too
   ;; abs(p) because p is allowed to be negative to eliminate bias
   if not (ss.planet[i].e.value lt (1d0-1d0/ss.planet[i].ar.value-abs(ss.planet[i].p.value)/ss.planet[i].ar.value)) then begin
      if keyword_set(verbose) then begin
         printandlog, 'Planet ' + strtrim(i,2) + ' will collide with the star!'+$
                      'e= ' + strtrim(ss.planet[i].e.value,2) + $
                      '; a/Rstar=' + strtrim(ss.planet[i].ar.value,2) + $
                      '; Rp/Rstar=' + strtrim(ss.planet[i].p.value,2) + $
                      '; Rstar=' + strtrim(ss.star.rstar.value,2), logname
         printandlog, "a/Rstar is derived using Kepler's law; adjust starting values " + $
                      'for mstar, mp, or period to change it', logname
      endif
      return, -1
   endif

   ;; error checking -- this should never happen
   if ~finite(ss.planet[i].ar.value) then begin
      printandlog, 'a/rstar is not finite; this should never happen', ss.logname
      stop
   endif

   if ss.planet[i].chord.fit then begin
      ;; derive cosi from chord
      ss.planet[i].b.value = sqrt((1d0+ss.planet[i].p.value)^2-ss.planet[i].chord.value^2)
      ss.planet[i].cosi.value = ss.planet[i].b.value/$
                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))

;***** debugging ********
      if ss.planet[i].cosi.value lt 0d0 then begin
         printandlog, string(ss.planet[i].b.value,format='("b=",f0.20)'), logname
         printandlog, string(ss.planet[i].ar.value,format='("a/r=",f0.20)'), logname
         printandlog, string(ss.planet[i].e.value,format='("e=",f0.20)'), logname
         printandlog, string(ss.planet[i].esinw.value,format='("esinw=",f0.20)'), logname
         printandlog, string(ss.planet[i].sign.value,format='("sign=",f0.20)'), logname
         printandlog, string(ss.planet[i].vcve.value,format='("vcve=",f0.20)'), logname
         printandlog, string(ss.planet[i].lsinw.value,format='("lsinw=",f0.20)'), logname
         printandlog, string(ss.planet[i].lsinw2.value,format='("lsinw2=",f0.20)'), logname
         printandlog, string(ss.planet[i].lcosw.value,format='("lcosw=",f0.20)'), logname
      endif
;*************************

   endif else if ss.planet[i].b.fit then begin
      ;; derive cosi from b
      ss.planet[i].cosi.value = ss.planet[i].b.value/$
                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))
   endif else if ss.planet[i].cosi.fit then begin
      ;; derive b from cosi
      ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*$
                             (1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; eq 7, Winn 2010
   endif
   ss.planet[i].bs.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)  ;; eq 8, Winn 2010

   ;; TC2TT can fail at weird inclinations; do this rejection here (instead of exofast_chi2v2.pro) 
   if ss.fittran[i] and (abs(ss.planet[i].b.value) gt (1d0 + ss.planet[i].p.value)) and (abs(ss.planet[i].bs.value) gt (1d0 + ss.planet[i].p.value)) then begin
      if keyword_set(verbose) then printandlog, 'planet does not transit; rejecting step', logname
      return, -1
   endif else if ss.planet[i].cosi.value lt 0 and ~ss.planet[i].i180 then begin
      if keyword_set(verbose) then printandlog, 'cosi out of bounds; rejecting step', logname
      return, -1
   endif else if abs(ss.planet[i].cosi.value) gt 1d0 then begin
      if keyword_set(verbose) then printandlog, 'cosi out of bounds; rejecting step', logname
      return, -1
   endif

   ;; derive inclination
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi

   ;; derive RV semi-amplitude
   ss.planet[i].k.value = (2d0*!dpi*G/(ss.planet[i].period.value*(ss.star.mstar.value + ss.planet[i].mpsun.value)^2d0))^(1d0/3d0) * $
                          ss.planet[i].mpsun.value*sin(ss.planet[i].i.value)/sqrt(1d0-ss.planet[i].e.value^2)*ss.constants.rsun/ss.constants.meter/ss.constants.day ;; m/s
   
   ;; Faigler & Mazeh, 2011
   ;; this assumes alpha_beam = 1d0, when in reality it is between 0.8
   ;; and 1.2
   if ss.planet[i].beam.derive and ~ss.planet[i].beam.fit then $
      ss.planet[i].beam.value = 4d0*ss.planet[i].k.value/(ss.constants.c/ss.constants.meter) ;; eq 1
;   if ss.planet[i].ellipsoidal.derive then $
;      ss.planet[i].ellipsoidal.value = ss.planet[i].mpsun.value/ss.star.mstar.value*sin(ss.planet[i].i.value)/ss.planet.ar.value^3*sin(ss.planet[i].i.value) ;; eq 2

   ;; derive lambda and bigomega
   ss.planet[i].lambda.value = atan(ss.planet[i].lsinlambda.value,ss.planet[i].lcoslambda.value)
   ss.planet[i].bigomega.value = atan(ss.planet[i].lsinbigomega.value,ss.planet[i].lcosbigomega.value)

   if ss.planet[i].tt.fit then begin

      ;; compute the time of conjunction
      ss.planet[i].tc.value = tc2tt(ss.planet[i].tt.value,$
                                    ss.planet[i].e.value,$
                                    ss.planet[i].i.value,$
                                    ss.planet[i].omega.value,$
                                    ss.planet[i].period.value,/reverse_correction)

   endif else begin
      ;; tc is fit, don't bother computing tt unless it's constrained by a prior
      if ss.planet[i].tt.prior ne 0 then begin
         
         ;; compute the time of minimum projected separation
         ss.planet[i].tt.value = tc2tt(ss.planet[i].tc.value,$
                                       ss.planet[i].e.value,$
                                       ss.planet[i].i.value,$
                                       ss.planet[i].omega.value,$
                                       ss.planet[i].period.value)
         
      endif
   endelse

   ;; time of periastron (required for the model)
   ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)
   ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*ss.planet[i].phase.value
   ;; pick the value closest to Tc
   if ss.planet[i].tp.value gt ss.planet[i].tc.value + ss.planet[i].period.value/2d0 then ss.planet[i].tp.value -= ss.planet[i].period.value
   if ss.planet[i].tp.value lt ss.planet[i].tc.value - ss.planet[i].period.value/2d0 then ss.planet[i].tp.value -= ss.planet[i].period.value

   ;; for plotting (fitting?) secondaries
   phase2 = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/secondary)
   ss.planet[i].ts.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phase2)

   ;; if given a prior, select the epoch closest to that prior
   if ss.planet[i].ts.prior ne 0d0 then begin
      nper = round((ss.planet[i].ts.prior - ss.planet[i].ts.value)/ss.planet[i].period.value)
      ss.planet[i].ts.value += nper*ss.planet[i].period.value
   endif else begin
      ;; otherwise select the epoch closest to Tc
      nper = round((ss.planet[i].tc.value - ss.planet[i].ts.value)/ss.planet[i].period.value)
      ss.planet[i].ts.value += nper*ss.planet[i].period.value
   endelse

endfor

for i=0L, ss.nband-1 do begin  
   massfraction = ss.planet[0].mpsun.value/(ss.star.mstar.value + ss.planet[0].mpsun.value)
   fluxfraction = ss.band[i].dilute.value
   ss.band[i].phottobary.value = 1d0/(massfraction-fluxfraction)
endfor

return, 1
end
