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

;; if value is a map, and the variable is fixed, we must propagate it here
for i=0L, n_elements(*ss.priors)-1 do begin

   ;; Skip priors with penalties -- the penalties will guide them to agree
   prior = (*ss.priors)[i]

   ;; if it's not fixed to a linked parameter, skip it. The penalties will sort it out
   if prior.gaussian_width ne 0d0 or prior.value[1] eq -1 then continue
   
   ;; the prior is linked to another variable -- get its value
   if prior.value[4] ne -1 then begin
      value = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].(prior.value[4])[prior.value[5]].value
   endif else if prior.value[2] ne -1 then begin
      value = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].value
   endif else begin
      value = ss.(prior.value[0])[prior.value[1]].value
   endelse

   ;; assign the value
   if prior.map[4] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]].(prior.map[4])[prior.map[5]].value = value
   endif else if prior.map[2] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]].value = value
   endif else if prior.map[0] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]].value = value
   endif     

endfor

;; derive stellar parameters
for i=0L, ss.nstars-1 do begin
   ss.star[i].mstar.value = 10^ss.star[i].logmstar.value
   ss.star[i].logg.value = alog10(ss.constants.gravitysun*ss.star[i].mstar.value/(ss.star[i].rstar.value^2)) ;; cgs
;; derive the distance from lstar
   ss.star[i].lstar.value = 4d0*!dpi*ss.star[i].rstar.value^2*ss.star[i].teff.value^4*sigmaB                                ;; L_sun
   if ss.star[i].distance.fit then ss.star[i].parallax.value = 1d3/ss.star[i].distance.value                                ;; mas
   if ss.star[i].parallax.fit then ss.star[i].distance.value = 1d3/ss.star[i].parallax.value                                ;; pc
   ss.star[i].fbol.value = (ss.star[i].lstar.value*ss.constants.lsun)/(4d0*!dpi*(ss.star[i].distance.value*ss.constants.pc)^2) ;; cgs
   ss.star[i].rhostar.value = ss.star[i].mstar.value/(ss.star[i].rstar.value^3)*ss.constants.rhosun                            ;; rho_sun
   ss.star[i].absks.value = ss.star[i].appks.value - 2.5d0*alog10((ss.star[i].distance.value/10d0)^2)                       ;; mag
endfor   

for j=0, ss.ntel-1 do begin
   if ss.telescope[j].jittervar.value gt 0 then $
      ss.telescope[j].jitter.value = sqrt(ss.telescope[j].jittervar.value) $
   else ss.telescope[j].jitter.value = 0d0
endfor

for i=0, ss.nplanets-1 do begin

   ;; derive the mass of the planet
   if ss.planet[i].logmp.fit then ss.planet[i].mpsun.value = 10^ss.planet[i].logmp.value $ ;; m_sun
   else ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value) ;; m_sun
   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup ;; m_jupiter
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth ;; m_earth

   ;; derive the radius of the planet   
   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star[ss.planet[i].starndx].rstar.value ;; r_sun
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup ;; r_jupiter
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth ;; r_earth

   ;; allow neg planet mass to avoid boundary bias, but can't have neg total mass
   if ss.star[ss.planet[i].starndx].mstar.value + ss.planet[i].mpsun.value le 0d0 then begin
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

      ;; what sign of the quadratic solution?
;      if ss.planet[i].sign2.fit then begin
      if ss.planet[i].sign.fit then begin
         ;; we fit for it
         sign = floor(ss.planet[i].sign.value)
      endif else begin
         ;; L implicitly defines it
         lsq = ss.planet[i].lsinw.value^2 + ss.planet[i].lcosw.value^2
         sign = (lsq ge 0.5d0) ;; if L^2 >= 0.5, use the positive solution to the quadratic!
      endelse

if 0 then begin
      ;; this ensures continuity across the vcve=1 boundary 
      ;; (as much as physically possible)
      if ss.planet[i].lsinw2.fit then begin
         if ss.planet[i].vcve.value gt 1d0 and ss.planet[i].lsinw2.value ge 0 and ~sign then begin
            ss.planet[i].lsinw.value = -ss.planet[i].lsinw2.value
            ss.planet[i].sign.value = ss.planet[i].sign2.value + 1d0
         endif else if ss.planet[i].vcve.value le 1d0 and ss.planet[i].lsinw2.value lt 0 and sign then begin
            ss.planet[i].lsinw.value = -ss.planet[i].lsinw2.value
            ss.planet[i].sign.value = ss.planet[i].sign2.value - 1d0
         endif else begin 
            ss.planet[i].lsinw.value = ss.planet[i].lsinw2.value
            ss.planet[i].sign.value = ss.planet[i].sign2.value
         endelse
      endif
endif else begin
;   ss.planet[i].lsinw.value = ss.planet[i].lsinw2.value
;   ss.planet[i].sign.value = ss.planet[i].sign2.value
endelse

      ;; based on L*cos(omega) and L*sin(omega)
;if (ss.planet[i].lsinw2.fit or ss.planet[i].lsinw.fit) and ss.planet[i].lcosw.fit then begin
      if ss.planet[i].lsinw.fit and ss.planet[i].lcosw.fit then begin
         ss.planet[i].omega.value = atan(ss.planet[i].lsinw.value, ss.planet[i].lcosw.value)
      endif

      ss.planet[i].e.value = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=floor(ss.planet[i].sign.value))

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
   ss.planet[i].arsun.value=(G*(ss.star[ss.planet[i].starndx].mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                       ;; semi-major axis in r_sun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star[ss.planet[i].starndx].rstar.value    ;; a/rstar (unitless)
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
                      '; Rstar=' + strtrim(ss.star[ss.planet[i].starndx].rstar.value,2), logname
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
if 0 then begin
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
   ;; if fitting TT, a primary transit is required (TT->TC conversion may fail if there is no transit). 
   ;; Otherwise, we need a primary or secondary transit (or the MCMC will be hopelessly lost)
   if ss.fittran[i] and (abs(ss.planet[i].b.value) gt (1d0 + ss.planet[i].p.value)) and ((abs(ss.planet[i].bs.value) gt (1d0 + ss.planet[i].p.value) or ss.planet[i].tt.fit)) then begin
      if keyword_set(verbose) then printandlog, 'planet ' + strtrim(i,2) + ' does not transit; rejecting step', logname
      return, -1
   endif else if ss.planet[i].cosi.value lt 0 and ~ss.planet[i].i180 then begin
      if keyword_set(verbose) then printandlog, 'cosi for planet ' + strtrim(i,2) + ' out of bounds; rejecting step', logname
      return, -1
   endif else if abs(ss.planet[i].cosi.value) gt 1d0 then begin
      if keyword_set(verbose) then printandlog, 'cosi for planet ' + strtrim(i,2) + ' out of bounds; rejecting step', logname
      return, -1
   endif

   ;; derive inclination
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi

   ;; derive RV semi-amplitude
   ss.planet[i].k.value = (2d0*!dpi*G/(ss.planet[i].period.value*(ss.star[ss.planet[i].starndx].mstar.value + ss.planet[i].mpsun.value)^2d0))^(1d0/3d0) * $
                          ss.planet[i].mpsun.value*sin(ss.planet[i].i.value)/sqrt(1d0-ss.planet[i].e.value^2)*ss.constants.rsun/ss.constants.meter/ss.constants.day ;; m/s
   
   ;; Faigler & Mazeh, 2011
   ;; this assumes alpha_beam = 1d0, when in reality it is between 0.8
   ;; and 1.2
   if ss.planet[i].beam.derive and ~ss.planet[i].beam.fit then $
      ss.planet[i].beam.value = 4d0*ss.planet[i].k.value/(ss.constants.c/ss.constants.meter) ;; eq 1
;   if ss.planet[i].ellipsoidal.derive then $
;      ss.planet[i].ellipsoidal.value = ss.planet[i].mpsun.value/ss.star[ss.planet[i].starndx].mstar.value*sin(ss.planet[i].i.value)/ss.planet.ar.value^3*sin(ss.planet[i].i.value) ;; eq 2

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

;for i=0L, ss.nband-1 do begin  
;   massfraction = ss.planet[0].mpsun.value/(ss.star[ss.band[i].starndx].mstar.value + ss.planet[0].mpsun.value)
;   fluxfraction = ss.band[i].dilute.value
;   ss.band[i].phottobary.value = 1d0/(massfraction-fluxfraction)
;endfor

return, 1
end
