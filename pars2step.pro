;; this function translates the starting guesses into the step parameters
function pars2step, ss

au = ss.constants.au/ss.constants.rsun
G = ss.constants.GMsun/ss.constants.RSun^3*ss.constants.day^2 ;; R_Sun^3/(m_sun*day^2)
mjup = ss.constants.gmjupiter/ss.constants.gmsun ;; m_sun
rjup = ss.constants.rjupiter/ss.constants.rsun  ;; r_sun
mearth = ss.constants.gmearth/ss.constants.gmsun ;; m_sun
rearth = ss.constants.rearth/ss.constants.rsun  ;; r_sun
sigmaB = ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; Stefan-Boltzmann constant

;; derive the stellar mass
if ss.star.mstar.userchanged then begin
   ss.star.logmstar.value = alog10( ss.star.mstar.value )
   if ss.star.mstar.scale ne 0d0 then $
      ss.star.logmstar.scale = ss.star.mstar.scale/(alog(10d0)*ss.star.mstar.value)
endif
ss.star.mstar.value = 10^ss.star.logmstar.value

;; derive teff, teffsed
if ss.star.teffsed.userchanged and ~ss.star.teff.userchanged then ss.star.teff.value = ss.star.teffsed.value
if ss.star.teff.userchanged and ~ss.star.teffsed.userchanged then ss.star.teffsed.value = ss.star.teff.value

;; derive rstar, rstarsed
if ss.star.rstarsed.userchanged and ~ss.star.rstar.userchanged then begin
   ss.star.rstar.value = ss.star.rstarsed.value
   ss.star.rstar.userchanged = 1B
endif
if ss.star.rstar.userchanged and ~ss.star.rstarsed.userchanged then ss.star.rstarsed.value = ss.star.rstar.value

;; derive the distance from the parallax, if given
if ~ss.star.distance.userchanged and ss.star.parallax.userchanged then $
   ss.star.distance.value = 1d3/ss.star.parallax.value

;; if a starting value for the stellar radius was given, 
;; use it to derive a starting point for the age
if ss.star.rstar.userchanged and ss.yy then begin
   ntries = 100
   age = dindgen(ntries)/(ntries-1)*13.82
   rstars = dblarr(ntries)
   for i=0, ntries-1 do begin
      junk = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, age[i], $
                            ss.star.teff.value,yyrstar=rstar, $
                            sigmab=ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2, $
                            gravitysun=ss.constants.gravitysun)
      if ~finite(junk) then stop;return, 0
      rstars[i] = rstar
   endfor
   junk = min(abs(rstars-ss.star.rstar.value),ndx)
   ss.star.age.value = age[ndx]
endif else if ss.yy then begin
   ;; otherwise derive the stellar radius
   junk = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, ss.star.age.value, $
                         ss.star.teff.value,yyrstar=rstar, $
                         sigmab=ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2, $
                         gravitysun=ss.constants.gravitysun)
   if ~finite(junk) then stop;return, 0
   ss.star.rstar.value = rstar
endif

;; if not specified, refine the starting value for eep (max out at 808)
if (ss.mist or ss.parsec) and ~ss.star.eep.userchanged and (ss.star.mstar.userchanged or ss.star.rstar.userchanged or ss.star.teff.userchanged) then begin  
   maxeep = 808d0
   mineep = 0d0
   neep = 30d0
   eeps = mineep + (maxeep-mineep)/(neep-1)*dindgen(neep)
   chi2 = dblarr(neep) + !values.d_infinity
   ages = dblarr(neep) + !values.d_infinity

   for i=0L, neep-1 do begin
      if ss.mist then begin
         mistchi2 =  massradius_mist(eeps[i],ss.star.mstar.value,ss.star.initfeh.value,$
                                     ss.star.age.value,ss.star.teff.value,$
                                     ss.star.rstar.value,ss.star.feh.value,mistage=mistage,/allowold)
         if finite(mistchi2) then begin
            if mistage gt 13.82d0 then break
            chi2[i] = mistchi2
            ages[i] = mistage
         endif
      endif else if ss.parsec then begin
         parsecchi2 =  massradius_parsec(eeps[i],ss.star.mstar.value,ss.star.initfeh.value,$
                                         ss.star.age.value,ss.star.teff.value,$
                                         ss.star.rstar.value,ss.star.feh.value,$
                                         parsec_age=parsec_age,/allowold)
         if finite(parsecchi2) then begin
            if parsec_age gt 13.82d0 then break
            chi2[i] = parsecchi2
            ages[i] = parsec_age
         endif
      endif
   endfor
   minchi2 = min(chi2,ndx)
   if ~finite(minchi2) then begin
      printandlog, 'No EEP between 0 and 808 produces a physical star. Adjust starting stellar parameters', ss.logname
      stop
   endif
   ss.star.eep.value = eeps[ndx]
   if ~ss.star.age.userchanged then ss.star.age.value = ages[ndx]

;; why did I do this? this was reallllly dumb
;; for high mass stars, the default age (4.6 Gyr) will select stars at
;; the end of their lifetimes (capped at carbon burning)
;   repeat begin
;      eep = (mineep + maxeep)/2d0
;      chi2 = massradius_mist(eep,ss.star.mstar.value,ss.star.initfeh.value,$
;                             ss.star.age.value,ss.star.teff.value,$
;                             ss.star.rstar.value,ss.star.feh.value,mistage=mistage,/allowold)
;      if ~finite(chi2) then begin
;         maxeep -= 1
;      endif else begin
;         if mistage gt ss.star.age.value then maxeep = eep $
;         else mineep = eep
;      endelse
;   endrep until abs(maxeep - mineep) lt 0.01
;   if ~finite(chi2) then stop; return, 0 
;   ss.star.eep.value = eep
endif
;; calculate the span of the data
minbjd=!values.d_infinity
maxbjd=-!values.d_infinity
if ss.ntran eq 0 then begin
   for i=0L, ss.ntel-1 do begin
      mint = min((*ss.telescope[i].rvptrs).bjd,max=maxt)
      if mint lt minbjd then minbjd = mint
      if maxt gt maxbjd then maxbjd = maxt
   endfor
endif else begin
   for i=0L, ss.ntran-1 do begin
      mint = min((*ss.transit[i].transitptrs).bjd,max=maxt)
      if mint lt minbjd then minbjd = mint
      if maxt gt maxbjd then maxbjd = maxt
   endfor
endelse

if ~ss.star.vsini.userchanged and (where(ss.doptom.dtscale.fit))[0] ne -1 then begin
   printandlog, 'You must set a prior on vsini', ss.logname
endif

;; derive quantities we'll use later
for i=0, ss.nplanets-1 do begin

   ;; derive the period from log(period)
   if ss.planet[i].period.userchanged then begin
      ss.planet[i].logp.value = alog10(ss.planet[i].period.value)
      
      ;; translate the stepping scale from P to logP
      if ss.planet[i].period.scale ne 0d0 then $
         ss.planet[i].logp.scale = ss.planet[i].period.scale/(alog(10d0)*ss.planet[i].period.value)
   endif else if ss.planet[i].logp.userchanged then begin      
      ss.planet[i].period.value = 10^ss.planet[i].logp.value
   endif else begin
      printandlog, 'You must set the starting value for the period', ss.logname
      stop
   endelse

   ;; how did the user seed e & omega?
   if ss.planet[i].secosw.userchanged or ss.planet[i].sesinw.userchanged then begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
      ss.planet[i].e.userchanged = 1B
      ss.planet[i].omega.userchanged = 1B
   endif else if ss.planet[i].qecosw.userchanged or ss.planet[i].qesinw.userchanged then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
      ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)      
      ss.planet[i].e.userchanged = 1B
      ss.planet[i].omega.userchanged = 1B
   endif else if ss.planet[i].ecosw.userchanged or ss.planet[i].esinw.userchanged then begin
      ss.planet[i].e.value = sqrt(ss.planet[i].ecosw.value^2 + ss.planet[i].esinw.value^2)
      ss.planet[i].omega.value = atan(ss.planet[i].esinw.value,ss.planet[i].ecosw.value)
      ss.planet[i].e.userchanged = 1B
      ss.planet[i].omega.userchanged = 1B
   endif



   ;; how did the user seed omega?
   if ss.planet[i].omega.userchanged then begin
      ;; do nothing
   endif else if ss.planet[i].omegadeg.userchanged then begin
      ss.planet[i].omega.value = ss.planet[i].omegadeg.value*!pi/180d0
      if ss.planet[i].omegadeg.scale ne 0d0 then $
         ss.planet[i].omega.scale = ss.planet[i].omegadeg.scale*!pi/180d0
   endif else if ss.planet[i].lsinw.userchanged or ss.planet[i].lsinw2.userchanged or ss.planet[i].lcosw.userchanged then begin
      if ss.planet[i].lsinw2.userchanged then begin
         if ss.planet[i].vcve.value lt 1d0 then ss.planet[i].lsinw.value = -ss.planet[i].lsinw2.value $
         else ss.planet[i].lsinw.value = ss.planet[i].lsinw2.value
      endif
      ss.planet[i].omega.value = atan(ss.planet[i].lsinw.value,ss.planet[i].lcosw.value)
      ss.planet[i].omega.userchanged = 1B
   endif

   ;; how did the user seed e?
   if ss.planet[i].e.userchanged then begin
      ;; do nothing
   endif else if ss.planet[i].vcve.userchanged then begin
      ;; pick the omega that maximizes the range of vcve
      if ~ss.planet[i].omega.userchanged then begin
         if ss.planet[i].vcve.value le 1d0 then ss.planet[i].omega.value = !dpi/2d0 $
         else ss.planet[i].omega.value = -!dpi/2d0
         ss.planet[i].omega.userchanged = 1B
      endif

      if ss.planet[i].sign.userchanged then begin
         sign = floor(ss.planet[i].sign.value)
         if ss.planet[i].vcve.value lt 1d0 then sign = ~sign
      endif

      ;; derive e (pick the sign if not chosen)
      ss.planet[i].e.value = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=sign)
      if ~ss.planet[i].sign.userchanged then begin
         sign = floor(sign)
         if ss.planet[i].vcve.value lt 1d0 then sign = ~sign
         ss.planet[i].sign.value = sign + 0.5d0
      endif

      ;; rederive lsinw/lcosw
      if ~ss.planet[i].lsinw.userchanged and ~ss.planet[i].lcosw.userchanged then begin
         if floor(ss.planet[i].sign.value) then L = sqrt(0.75d0) $
         else L = sqrt(0.25d0)
         ss.planet[i].lsinw.value = L*sin(ss.planet[i].omega.value)
         if ss.planet[i].vcve.value - 1d0 lt 0d0 then begin
            ss.planet[i].lsinw2.value = -L*sin(ss.planet[i].omega.value)
         endif else ss.planet[i].lsinw2.value = L*sin(ss.planet[i].omega.value)
         ss.planet[i].lcosw.value = L*cos(ss.planet[i].omega.value)
      endif

   endif 

   ;; error checking on e/omega
   if ss.planet[i].e.value ge 1d0 or ss.planet[i].e.value lt 0d0 or ~finite(ss.planet[i].e.value) then stop;return, 0
   if ~finite(ss.planet[i].omega.value) then stop;return, 0

   ;; translate from e/omega to various e/omega parameterizations
   if ~ss.planet[i].qecosw.userchanged then ss.planet[i].qecosw.value = (ss.planet[i].e.value)^(0.25d0)*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].qesinw.userchanged then ss.planet[i].qesinw.value = (ss.planet[i].e.value)^(0.25d0)*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].secosw.userchanged then ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].sesinw.userchanged then ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].ecosw.userchanged then ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].esinw.userchanged then ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].vcve.userchanged then ss.planet[i].vcve.value = sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; need to set the sign to recover desired e
   if ~ss.planet[i].sign.userchanged then begin

      e1 = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=0B)
      e2 = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=1B)

      if abs(ss.planet[i].e.value - e1) lt 1d-8 then ss.planet[i].sign.value = 0.5 $
      else if abs(ss.planet[i].e.value - e2) lt 1d-8 then ss.planet[i].sign.value = 1.5 $
      else begin
         printandlog, 'ERROR translating staring e and omega to parameterization', ss.logname
         stop
      endelse
      if ss.planet[i].vcve.value lt 1d0 then ss.planet[i].sign.value = ~floor(ss.planet[i].sign.value) + 0.5

   endif

   ;; rederive lsinw/lcosw
   if ~ss.planet[i].lsinw.userchanged and ~ss.planet[i].lcosw.userchanged then begin
      if floor(ss.planet[i].sign.value) then L = sqrt(0.75d0) $
      else L = sqrt(0.25d0)
      ss.planet[i].lsinw.value = L*sin(ss.planet[i].omega.value)
      if ss.planet[i].vcve.value -1d0 lt 0d0 then begin
         ss.planet[i].lsinw2.value = -L*sin(ss.planet[i].omega.value)
      endif else ss.planet[i].lsinw2.value = L*sin(ss.planet[i].omega.value)
      ss.planet[i].lcosw.value = L*cos(ss.planet[i].omega.value)
   endif

   ;; how did the user seed rp/rstar?
   if ss.planet[i].p.userchanged then begin
      if ~ss.planet[i].rpearth.userchanged then $
         ss.planet[i].rpearth.value = ss.planet[i].p.value*ss.star.rstar.value/rearth
      ss.planet[i].rpearth.userchanged = 1B
   endif else if ss.planet[i].rpearth.userchanged then begin
      ss.planet[i].p.value = ss.planet[i].rpearth.value*rearth/ss.star.rstar.value
   endif else if ss.planet[i].rp.userchanged then begin
      ss.planet[i].p.value = ss.planet[i].rp.value*rjup/ss.star.rstar.value
      if ~ss.planet[i].rpearth.userchanged then $
         ss.planet[i].rpearth.value = ss.planet[i].p.value*ss.star.rstar.value/rearth
      ss.planet[i].rpearth.userchanged = 1B
   endif else if ss.planet[i].mpearth.userchanged then begin 
      ss.planet[i].rpearth.value = massradius_chen(ss.planet[i].mpearth.value)
      ss.planet[i].p.value = ss.planet[i].rpearth.value*rearth/ss.star.rstar.value
   endif

   ;; how did the user seed the inclination?
   if ss.planet[i].i.userchanged then begin
      sini = sin(ss.planet[i].i.value)
   endif else if ss.planet[i].ideg.userchanged then begin      
     ss.planet[i].i.value = ss.planet[i].ideg.value*!dpi/180d0     
     sini = sin(ss.planet[i].i.value)
   endif else if ss.planet[i].chord.userchanged then begin
      if ~ss.planet[i].chord.fit then $
         printandlog, 'changing starting value for chord when not fitting chord is not supported', ss.logname
;      ;; need to know rp/rstar, a/rstar, mpsun, inc (!) first
;      ss.planet[i].b.value = sqrt((1d0+ss.planet[i].p.value)^2-ss.planet[i].chord.value^2)
;      ss.planet[i].cosi.value = ss.planet[i].b.value/$
;                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))
;      ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
      sini = 1d0
   endif else if ss.planet[i].b.userchanged then begin
      if ~ss.planet[i].b.fit then $
         printandlog, 'changing starting value for b when not fitting b is not supported', ss.logname
;      ;; need to know a/rstar, mpsun, inc (!) first
;      ss.planet[i].cosi.value = ss.planet[i].b.value/$
;                                (ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value))
;      ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
      sini = 1d0
   endif else begin
      ;; use the default cosi starting value
      ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
      sini = sin(ss.planet[i].i.value)
   endelse      
  
   ;; how did the user seed the planet mass?
   if ss.planet[i].mpsun.userchanged then begin
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
      ;; translate the stepping scale from mpsun to logmp
      if ss.planet[i].logmp.scale ne 0d0 then $
         ss.planet[i].logmp.scale = ss.planet[i].mpsun.scale/(alog(10d0)*ss.planet[i].mpsun.value)
   endif else if ss.planet[i].logmp.userchanged then begin
      ss.planet[i].mpsun.value = 10d0^ss.planet[i].logmp.value
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].mp.userchanged then begin
      ;; convert from Jupiter masses to Solar masses
      ss.planet[i].mpsun.value = ss.planet[i].mp.value*mjup
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].mpearth.userchanged then begin
      ;; convert from Earth masses to Solar masses
      ss.planet[i].mpsun.value = ss.planet[i].mp.value*mearth
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].msiniearth.userchanged then begin
      ;; convert from Earth masses to Solar masses
      ss.planet[i].mpsun.value = ss.planet[i].msiniearth.value*mearth/sini
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].msini.userchanged then begin
      ;; convert from Jupiter masses to Solar masses
      ss.planet[i].mpsun.value = ss.planet[i].msini.value*mjup/sini
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].logk.userchanged then begin
      ss.planet[i].k.value = 10d0^ss.planet[i].logk.value
      ;; now convert from K to mpsun
      ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                       ss.planet[i].i.value, ss.planet[i].period.value, $
                                       ss.star.mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].k.userchanged then begin
      if ss.planet[i].k.value lt 0d0 then begin
         ss.planet[i].mpsun.value = -ktom2(-ss.planet[i].K.value, ss.planet[i].e.value,$
                                           ss.planet[i].i.value, ss.planet[i].period.value, $
                                           ss.star.mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      endif else begin
         ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                          ss.planet[i].i.value, ss.planet[i].period.value, $
                                          ss.star.mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      endelse
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].rpearth.userchanged then begin
      ss.planet[i].mpearth.value = massradius_chenreverse(ss.planet[i].rpearth.value)
      ss.planet[i].mpsun.value = ss.planet[i].mpearth.value*mearth ;; m_sun
   endif

   ;; how did the user seed a/rstar?
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                       ;; semi-major axis in r_sun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value    ;; a/rstar (unitless)
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; semi major axis in AU
   ss.planet[i].cosi.scale = 1d0/ss.planet[i].ar.value

   ;; translate from bigomegadeg to bigomega
   if ss.planet[i].bigomegadeg.userchanged then begin
      ss.planet[i].bigomega.value = ss.planet[i].bigomegadeg.value*!pi/180d0
      if ss.planet[i].bigomegadeg.scale ne 0 then $
         ss.planet[i].bigomega.scale = ss.planet[i].bigomegadeg.scale*!pi/180d0
   endif

   ;; translate from bigomega to lsinbigomega/lcosbigomega
   if ss.planet[i].bigomega.userchanged then begin
      ss.planet[i].lsinbigomega.value = 0.5d0*sin(ss.planet[i].bigomega.value) 
      ss.planet[i].lcosbigomega.value = 0.5d0*cos(ss.planet[i].bigomega.value) 
   endif

   ;; translate from lambdadeg to lambda
   if ss.planet[i].lambdadeg.userchanged then begin
      ss.planet[i].lambda.value = ss.planet[i].lambdadeg.value*!pi/180d0
      if ss.planet[i].lambdadeg.scale ne 0 then $
         ss.planet[i].lambda.scale = ss.planet[i].lambdadeg.scale*!pi/180d0
   endif

   ;; translate from lambda to lsinlambda/lcoslambda
   if ss.planet[i].lambda.userchanged then begin
      ss.planet[i].lsinlambda.value = sin(ss.planet[i].lambda.value) 
      ss.planet[i].lcoslambda.value = cos(ss.planet[i].lambda.value) 
   endif

   ;; translate to tt
   if ss.planet[i].tt.fit then begin
      if ss.planet[i].tt.userchanged then begin
         ;; do nothing; it was already set explicitly via the prior file
      endif else if ss.planet[i].tc.userchanged then begin
         ;; derive tt from tc
         ss.planet[i].tt.value = tc2tt(ss.planet[i].tc.value,$
                                       ss.planet[i].e.value,$
                                       ss.planet[i].i.value,$
                                       ss.planet[i].omega.value,$
                                       ss.planet[i].period.value)

         ;; make tt as close as possible to tc
         nper = round((ss.planet[i].tc.value - ss.planet[i].tt.value)/ss.planet[i].period.value)
         ss.planet[i].tt.value += nper*ss.planet[i].period.value

      endif else if ss.planet[i].tp.userchanged then begin
         ;; derive tp->tc->tt
         phase = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value, /primary)
         ss.planet[i].tc.value = ss.planet[i].tp.value + ss.planet[i].period.value*phase
         ss.planet[i].tt.value = tc2tt(ss.planet[i].tc.value,$
                                       ss.planet[i].e.value,$
                                       ss.planet[i].i.value,$
                                       ss.planet[i].omega.value,$
                                       ss.planet[i].period.value)

         ;; make tt as close as possible to tp
         nper = round((ss.planet[i].tp.value - ss.planet[i].tt.value)/ss.planet[i].period.value)
         ss.planet[i].tt.value += nper*ss.planet[i].period.value

         ;; make tc as close as possible to tp
         nper = round((ss.planet[i].tp.value - ss.planet[i].tc.value)/ss.planet[i].period.value)
         ss.planet[i].tc.value += nper*ss.planet[i].period.value

      endif else begin
         ;; otherwise guess it's in the middle of the data
         ss.planet[i].tt.value = (maxbjd + minbjd)/2d0
         printandlog, 'Warning: No starting guess for TT (or TC or TP) given; assuming the middle of the data', ss.logname
      endelse
      if ss.planet[i].tt.prior eq 0d0 then ss.planet[i].tt.prior = ss.planet[i].tt.value 
   endif else begin
      ;; fit tc
      if ss.planet[i].tc.userchanged then begin
         ;; do nothing; it was already set explicitly via the prior file
      endif else if ss.planet[i].tp.userchanged then begin
         ;; derive tp->tc
         phase = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value, /primary)
         ss.planet[i].tc.value = ss.planet[i].tp.value + ss.planet[i].period.value*phase
         ss.planet[i].tc.prior = ss.planet[i].tc.value 
      endif else if ss.planet[i].tt.userchanged then begin
         ss.planet[i].tc.value = tc2tt(ss.planet[i].tt.value,$
                                       ss.planet[i].e.value,$
                                       ss.planet[i].i.value,$
                                       ss.planet[i].omega.value,$
                                       ss.planet[i].period.value,/reverse_correction)

         ;; make tt as close as possible to tc
         nper = round((ss.planet[i].tt.value - ss.planet[i].tc.value)/ss.planet[i].period.value)
         ss.planet[i].tc.value += nper*ss.planet[i].period.value
      endif else begin
         ;; otherwise guess it's in the middle of the data
         ss.planet[i].tc.value = (maxbjd + minbjd)/2d0
         printandlog, 'Warning: No starting guess for TC (or TT or TP) given; assuming the middle of the data', ss.logname
      endelse
      if ss.planet[i].tc.prior eq 0d0 then ss.planet[i].tc.prior = ss.planet[i].tc.value 
   endelse

   ;; approximate duration taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)  ;; eq 7, Winn 2010
   ss.planet[i].t14.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+abs(ss.planet[i].p.value))^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)


   if ~ss.planet[i].chord.userchanged then ss.planet[i].chord.value = sqrt((1d0+ss.planet[i].p.value)^2 - ss.planet[i].b.value^2)

endfor

return, 1

end
