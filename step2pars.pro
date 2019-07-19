function step2pars, ss, verbose=verbose, logname=logname, changedefaults=changedefaults

G = ss.constants.GMSun/ss.constants.RSun^3*ss.constants.day^2 ;; R_sun^3/(m_sun*day^2)
AU = ss.constants.au/ss.constants.rsun ;; R_sun
mjup = ss.constants.gmjupiter/ss.constants.gmsun ;; m_sun
rjup = ss.constants.rjupiter/ss.constants.rsun  ;; r_sun
mearth = ss.constants.gmearth/ss.constants.gmsun ;; m_sun
rearth = ss.constants.rearth/ss.constants.rsun  ;; r_sun
sigmaB = ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; Stefan-Boltzmann constant

mindist = dblarr(ss.nplanets>1)
maxdist = dblarr(ss.nplanets>1)

;; derive stellar parameters
ss.star.mstar.value = 10^ss.star.logmstar.value
ss.star.logg.value = alog10(ss.constants.gravitysun*ss.star.mstar.value/(ss.star.rstar.value^2)) ;; cgs
;; derive the distance from lstar
ss.star.lstar.value = 4d0*!dpi*ss.star.rstar.value^2*ss.star.teff.value^4*sigmaB ;; L_sun

if ss.star.distance.fit then ss.star.parallax.value = 1d3/ss.star.distance.value ;; mas
if ss.star.parallax.fit then ss.star.distance.value = 1d3/ss.star.parallax.value ;; mas
;ss.star.fbol.value = (ss.star.lstar.value/ss.constants.lsun)/(4d0*!dpi*(ss.star.distance.value/ss.constants.pc)^2) ;; cgs

;print, 'derived: ' + strtrim(ss.star.fbol.value,2)
;stop

;; require planets to stay in the same order as they start
if ss.nplanets gt 0 then begin
   if max(abs(ss.planetorder - sort(ss.planet.logp.value))) ne 0 then begin
;      if keyword_set(verbose) then printandlog, 'Planets must remain in the original order! Rejecting step', logname
;      stop
;      return, -1
   endif
endif

for j=0, ss.ntel-1 do begin
   if ss.telescope[j].jittervar.value gt 0 then $
      ss.telescope[j].jitter.value = sqrt(ss.telescope[j].jittervar.value)
endfor

for i=0, ss.nplanets-1 do begin

   ;; derive K if logK is fit
   if ss.planet[i].logk.fit then ss.planet[i].k.value = 10^ss.planet[i].logk.value
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].period.value = 10^ss.planet[i].logp.value

   ;; bound tc to be within 1 period of the starting value for tc
   if ss.planet[i].tc.prior ne 0d0 then begin
      nper = round((ss.planet[i].tc.prior - ss.planet[i].tc.value)/ss.planet[i].period.value)
      ss.planet[i].tc.value -= nper*ss.planet[i].period.value
   endif

   ;; eccentricity and argument of periastron
   if ss.planet[i].secosw.fit and ss.planet[i].sesinw.fit then begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endif else if ss.planet[i].ecosw.fit and ss.planet[i].esinw.fit then begin
      ss.planet[i].e.value = sqrt(ss.planet[i].ecosw.value^2 + ss.planet[i].esinw.value^2)
      ss.planet[i].omega.value = atan(ss.planet[i].esinw.value,ss.planet[i].ecosw.value)
   endif else if ss.planet[i].e.fit and ss.planet[i].omega.fit then begin
      if ss.planet[i].e.value lt 0d0 then begin
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

   ;; bound lambda to be between -pi and pi
   ss.planet[i].lambda.value = ss.planet[i].lambda.value mod (2d0*!dpi)
   if ss.planet[i].lambda.value gt !dpi then ss.planet[i].lambda.value -= (2d0*!dpi)
   if ss.planet[i].lambda.value lt -!dpi then ss.planet[i].lambda.value += (2d0*!dpi)

   if ss.fitrv[i] or ss.fittran[i] or i gt 1 then begin
      ;; fitting transits, RVs, or multiple planets mean we can resolve the discreet degeneracy
      ;; between Omega and Omega + pi
      if ss.planet[i].bigomega.value ge (2d0*!dpi) then ss.planet[i].bigomega.value -= 2d0*!dpi
      if ss.planet[i].bigomega.value lt 0d0  then ss.planet[i].bigomega.value += 2d0*!dpi
   endif else begin
      ;; otherwise, we can only measure the longitude of the Node 
      ;; (0 < Omega < 180)
      if ss.planet[i].bigomega.value ge (!dpi) then ss.planet[i].bigomega.value -= !dpi
      if ss.planet[i].bigomega.value lt 0d0  then ss.planet[i].bigomega.value += !dpi
   endelse

   if ss.planet[i].k.value le 0d0 then begin
      ss.planet[i].mpsun.value = -ktom2(-ss.planet[i].K.value, ss.planet[i].e.value,$
                                        ss.planet[i].i.value, ss.planet[i].period.value, $
                                        ss.star.mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun      
   endif else ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                               ss.planet[i].i.value, ss.planet[i].period.value, $
                                               ss.star.mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun

;   if ss.planet[i].mpsun.value gt 0.08d0 then begin
;      if keyword_set(verbose) then printandlog, 'Planet ' + strtrim(i,2) + ' mass (' + strtrim(ss.planet[i].mpsun.value/mjup,2) + ' M_J) above hydrogen burning limit',logname
;      return, -1 ;; planet above the hydrogen burning limit
;   endif

   if ss.star.mstar.value + ss.planet[i].mpsun.value le 0d0 then begin
      if keyword_set(verbose) then printandlog, 'Net mass (Mstar + Mp) must be positive; rejecting step', logname
      return, -1
   endif

   if ~finite(ss.planet[i].mpsun.value) then begin
      if keyword_set(verbose) then printandlog, 'Planet mass must be finite; rejecting step', logname
      return, -1
   endif

   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup ;; m_jupiter
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth ;; m_earth
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                    ;; semi-major axis in r_sun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value ;; a/rstar (unitless)
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; semi major axis in AU

   ;; tighter (empirical) constraint on eccentricity (see Eastman 2013)
   if ss.tides and ss.planet[i].e.value gt (1d0-3d0/ss.planet[i].ar.value) then begin
      if keyword_set(verbose) then printandlog, 'TIDES flag set and eccentricity (' + strtrim(ss.planet[i].e.value,2) + ') exceeds empirical limit (' + strtrim((1d0-3d0/ss.planet[i].ar.value),2) + ')', logname
      return, -1
   endif ;else if ss.planet[i].e.value gt (1d0-1d0/ss.planet[i].ar.value)

   ;; limit eccentricity to avoid collision with star during periastron
   ;; the ignored tidal effects would become important long before this,
   ;; but this prevents numerical problems compared to the e < 1 constraint
   ;; the not/lt (instead of ge) robustly handles NaNs too
   ;; abs(p) because p is allowed to be negative to eliminate bias
   if not (ss.planet[i].e.value lt (1d0-1d0/ss.planet[i].ar.value-abs(ss.planet[i].p.value)/ss.planet[i].ar.value)) then begin
      if keyword_set(verbose) then printandlog, 'Planet ' + strtrim(i,2) + ' will collide with the star! e=' + strtrim(ss.planet[i].e.value,2) + '; a/Rstar=' + strtrim(ss.planet[i].ar.value,2) + '; Rp/Rstar=' + strtrim(ss.planet[i].p.value,2) + '; Rstar=' + strtrim(ss.star.rstar.value,2), logname
      if keyword_set(verbose) then printandlog, 'Rstar is derived from YY isochrones; adjust starting values for Age, Teff, [Fe/H], or Mstar to change it', logname
      return, -1
   endif

   if ss.planet[i].e.value eq 0d0 then ss.planet[i].omega.value = !dpi/2d0

   if ~finite(ss.planet[i].ar.value ) then stop

   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star.rstar.value ;; r_sun
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup ;; r_jupiter
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth ;; r_earth

   ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; eq 7, Winn 2010

   ;; no transit and we're fitting a transit! This causes major problems; exclude this model
   if ss.planet[i].fittran and (ss.planet[i].b.value gt (1d0+ss.planet[i].p.value)) then begin
       if keyword_set(verbose) then printandlog, 'Planet does not transit!', logname
      return, -1
   endif

   ;; time of periastron
   ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)
   ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*ss.planet[i].phase.value

   if keyword_set(changedefaults) then begin
      ;; if the user supplied a value rp/rstar but not K, use
      ;; the mass-radius relation to derive a better starting value for K
      if ss.planet[i].k.value eq 10d0 and ss.planet[i].p.value ne 0.1d0 then begin
         ss.planet[i].mpearth.value = massradius_chenreverse(ss.planet[i].rpearth.value)
         ss.planet[i].mpsun.value = ss.planet[i].mpearth.value*mearth ;; m_sun
         ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup        ;; m_jupiter
         ss.planet[i].k.value = (2d0*!dpi*G/(ss.planet[i].period.value*(ss.star.mstar.value + ss.planet[i].mpsun.value)^2d0))^(1d0/3d0) * $
                                ss.planet[i].mpsun.value*sin(ss.planet[i].i.value)/sqrt(1d0-ss.planet[i].e.value^2)*ss.constants.rsun/ss.constants.meter/ss.constants.day ;; m/s
         ss.planet[i].logk.value = alog10(ss.planet[i].k.value)
      endif

      ;; if the user supplied a value of K changed, but not rp/rstar, use
      ;; the mass-radius relation to derive a better starting value for rp/rstar
      if ss.planet[i].p.value eq 0.1d0 and ss.planet[i].k.value ne 10d0 then begin
         ss.planet[i].rpearth.value = massradius_chen(ss.planet[i].mpearth.value)
         ss.planet[i].rpsun.value = ss.planet[i].rpearth.value*rearth ;; r_sun
         ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup        ;; r_jupiter
         ss.planet[i].p.value = ss.planet[i].rpsun.value/ss.star.rstar.value  
      endif
   endif

endfor

for i=0L, ss.nband-1 do begin  
   massfraction = ss.planet[0].mpsun.value/(ss.star.mstar.value + ss.planet[0].mpsun.value)
   fluxfraction = ss.band[i].dilute.value
   ss.band[i].phottobary.value = 1d0/(massfraction-fluxfraction)
endfor

if ~ss.alloworbitcrossing then begin
   for i=0L, ss.nplanets-1 do begin
      ;; for multi-planet systems, make sure they don't enter each other's hill spheres
      ;; if mp unknown, mp=0 => hill radius=0 => planets can't cross orbits
      ;; **** ignores mutual inclination; a priori excludes systems like Neptune and Pluto!! ****
      hillradius = ((1d0-ss.planet[i].e.value)*ss.planet[i].a.value*(ss.planet[i].mpsun.value/(3d0*ss.star.mstar.value))^(1d0/3d0)) > 0d0
      mindist[i] = (1d0-ss.planet[i].e.value)*ss.planet[i].a.value - hillradius
      maxdist[i] = (1d0+ss.planet[i].e.value)*ss.planet[i].a.value + hillradius
      for j=i-1,0,-1 do begin
         if ((mindist[i] ge mindist[j]) and (mindist[i] le maxdist[j])) or $
            ((maxdist[i] ge mindist[j]) and (maxdist[i] le maxdist[j])) or $
            ((mindist[j] ge mindist[i]) and (mindist[j] le maxdist[i])) or $
            ((maxdist[j] ge mindist[i]) and (maxdist[j] le maxdist[i])) then begin
            if keyword_set(verbose) then printandlog, 'Planets ' + strtrim(j,2) + ' (' + strtrim(mindist[j],2) + ',' + strtrim(maxdist[j],2) + ') and ' + strtrim(i,2) + ' (' + strtrim(mindist[i],2) + ',' + strtrim(maxdist[i],2) + ') cross paths', logname
            return, -1
         endif
      endfor
   endfor
endif

return, 1
end
