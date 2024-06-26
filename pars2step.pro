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
for i=0, ss.nstars-1 do begin
   if ss.star[i].mstar.userchanged then begin
      ss.star[i].logmstar.value = alog10( ss.star[i].mstar.value )
      if ss.star[i].mstar.scale ne 0d0 then $
         ss.star[i].logmstar.scale = ss.star[i].mstar.scale/(alog(10d0)*ss.star[i].mstar.value)
   endif
   ss.star[i].mstar.value = 10^ss.star[i].logmstar.value

   ;; derive teff, teffsed, feh, fehsed
   if ss.star[i].teffsed.userchanged and ~ss.star[i].teff.userchanged then ss.star[i].teff.value = ss.star[i].teffsed.value
   if ss.star[i].teff.userchanged and ~ss.star[i].teffsed.userchanged then ss.star[i].teffsed.value = ss.star[i].teff.value
   if ss.star[i].fehsed.userchanged and ~ss.star[i].feh.userchanged then ss.star[i].feh.value = ss.star[i].fehsed.value
   if ss.star[i].feh.userchanged and ~ss.star[i].fehsed.userchanged then ss.star[i].fehsed.value = ss.star[i].feh.value
   if ss.star[i].feh.userchanged and ~ss.star[i].initfeh.userchanged then ss.star[i].initfeh.value = ss.star[i].feh.value
   
   ;; derive rstar, rstarsed
   if ss.star[i].rstarsed.userchanged and ~ss.star[i].rstar.userchanged then begin
      ss.star[i].rstar.value = ss.star[i].rstarsed.value
      ss.star[i].rstar.userchanged = 1B
   endif
   if ss.star[i].rstar.userchanged and ~ss.star[i].rstarsed.userchanged then ss.star[i].rstarsed.value = ss.star[i].rstar.value

   ;; derive the distance from the parallax, if given
   if ~ss.star[i].distance.userchanged and ss.star[i].parallax.userchanged then $
      ss.star[i].distance.value = 1d3/ss.star[i].parallax.value
   
   ;; if a starting value for the stellar radius was given, 
   ;; use it to derive a starting point for the age
   if ss.star[i].rstar.userchanged and ss.yy[i] then begin
      ntries = 100
      age = dindgen(ntries)/(ntries-1)*13.82
      rstars = dblarr(ntries)
      for j=0L, ntries-1 do begin
         junk = massradius_yy3(ss.star[i].mstar.value, ss.star[i].feh.value, age[i], $
                               ss.star[i].teff.value,yyrstar=rstar, $
                               sigmab=ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2, $
                               gravitysun=ss.constants.gravitysun)
         if ~finite(junk) then return, 0
         rstars[j] = rstar
      endfor
      junk = min(abs(rstars-ss.star[i].rstar.value),ndx)
      ss.star[i].age.value = age[ndx]
   endif else if ss.yy[i] then begin
      ;; otherwise derive the stellar radius
      junk = massradius_yy3(ss.star[i].mstar.value, ss.star[i].feh.value, ss.star[i].age.value, $
                            ss.star[i].teff.value,yyrstar=rstar, $
                            sigmab=ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2, $
                            gravitysun=ss.constants.gravitysun)
      if ~finite(junk) then return, 0
      ss.star[i].rstar.value = rstar
   endif
   
   ;; if not specified, refine the starting value for eep (max out at 808)
   if (ss.mist[i] or ss.parsec[i]) and ~ss.star[i].eep.userchanged and (ss.star[i].mstar.userchanged or ss.star[i].rstar.userchanged or ss.star[i].teff.userchanged) then begin  
      maxeep = 808d0
      mineep = 0d0
      neep = 100d0
      eeps = mineep + (maxeep-mineep)/(neep-1)*dindgen(neep)
      chi2 = dblarr(neep) + !values.d_infinity
      ages = dblarr(neep) + !values.d_infinity

      for j=0L, neep-1 do begin
         if ss.mist[i] then begin
            mistchi2 =  massradius_mist(eeps[j],ss.star[i].mstar.value,ss.star[i].initfeh.value,$
                                        ss.star[i].age.value,ss.star[i].teff.value,$
                                        ss.star[i].rstar.value,ss.star[i].feh.value,mistage=mistage,fitage=1,$
                                        /allowold,tefffloor=ss.teffemfloor,fehfloor=ss.fehemfloor,$
                                        rstarfloor=ss.rstaremfloor, agefloor=ss.ageemfloor)
            if eeps[j] lt 202 then mistchi2+=30 ;; penalize PMS stars
            if finite(mistchi2) then begin
               if mistage gt 13.82d0 then break
               chi2[j] = mistchi2
               ages[j] = mistage
            endif
         endif else if ss.parsec[i] then begin
            parsecchi2 =  massradius_parsec(eeps[j],ss.star[i].mstar.value,ss.star[i].initfeh.value,$
                                            ss.star[i].age.value,ss.star[i].teff.value,$
                                            ss.star[i].rstar.value,ss.star[i].feh.value,$
                                            parsec_age=parsec_age,fitage=1,/allowold,$
                                            tefffloor=ss.teffemfloor,fehfloor=ss.fehemfloor,$
                                            rstarfloor=ss.rstaremfloor, agefloor=ss.ageemfloor)
            if eeps[j] lt 202 then parsecchi2+=30 ;; penalize PMS stars
            if finite(parsecchi2) then begin
               if parsec_age gt 13.82d0 then break
               chi2[j] = parsecchi2
               ages[j] = parsec_age
            endif
         endif
      endfor
      minchi2 = min(chi2,ndx)
      if ~finite(minchi2) then begin
         printandlog, 'No EEP between 0 and 808 produces a physical star. Adjust starting stellar parameters', ss.logname
         stop
      endif
      ss.star[i].eep.value = eeps[ndx]
      if ~ss.star[i].age.userchanged then ss.star[i].age.value = ages[ndx]
      
   endif
endfor

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

if ~ss.star[0].vsini.userchanged and (where(ss.doptom.dtscale.fit))[0] ne -1 then begin
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
      ;; translate degrees to radians
      ss.planet[i].omega.value = ss.planet[i].omegadeg.value*!pi/180d0
      if ss.planet[i].omegadeg.scale ne 0d0 then $
         ss.planet[i].omega.scale = ss.planet[i].omegadeg.scale*!pi/180d0
   endif else if ss.planet[i].lsinw.userchanged or ss.planet[i].lcosw.userchanged then begin
      ;; translate Lcosw/Lsinw to omega
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

;      if ss.planet[i].sign2.userchanged then begin
      if ss.planet[i].sign.userchanged then begin
         sign = floor(ss.planet[i].sign.value)
      endif else undefine, sign ;; can't have it defined from previous planets or it will fail if it changes

      ;; rederive lsinw/lcosw
      if ~ss.planet[i].lsinw.userchanged and ~ss.planet[i].lcosw.userchanged then begin
         if floor(ss.planet[i].sign.value) then L = sqrt(0.75d0) $
         else L = sqrt(0.25d0)
         ss.planet[i].lsinw.value = L*sin(ss.planet[i].omega.value)
         ss.planet[i].lcosw.value = L*cos(ss.planet[i].omega.value)
      endif

      ;; derive e (pick the sign if not chosen)
      ss.planet[i].e.value = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=sign)
      ss.planet[i].sign.value = sign
      ss.planet[i].omegadeg.value = ss.planet[i].omega.value*180d0/!dpi

   endif else if ss.planet[i].tc.userchanged and ss.planet[i].ts.userchanged then begin
      ;; from eclipse timing
      print, 'warning: setting e/omega from times and durations is untested'
      if ss.planet[i].esinw.userchanged then esinw=ss.planet[i].esinw.value
      if ss.planet[i].tfwhm.userchanged then tfwhmp=ss.planet[i].tfwhm.value
      if ss.planet[i].tfwhms.userchanged then tfwhms=ss.planet[i].tfwhms.value
      if ss.planet[i].i.userchanged then inc=ss.planet[i].i.value
      if ss.planet[i].p.userchanged then p=ss.planet[i].p.value
      ;; this is approximate. all other ways to initialize should take precedence
      junk = time2ew(ss.planet[i].tc.value, ss.planet[i].ts.value,$
                     ss.planet[i].period.value,esinwin=esinw,$
                     tfwhmp=tfwhmp, tfwhms=tfwhms, p=p, inc=inc,$
                     ecc=e,omega=omega,ecosw=ecosw)
      if ~ss.planet[i].e.userchanged then ss.planet[i].e.value = e
      if ~ss.planet[i].omega.userchanged then ss.planet[i].omega.value = omega
      ss.planet[i].e.userchanged = 1B
      ss.planet[i].omega.userchanged = 1B
   endif 

   ;; error checking on e/omega
   if ss.planet[i].e.value ge 1d0 or ss.planet[i].e.value lt 0d0 or ~finite(ss.planet[i].e.value) then begin
      printandlog, 'The starting value for eccentricity_' + strtrim(i,2) + ' is not physical. If you are seeing this with an automatically generated prior file after updating EXOFASTv2, delete the starting points for sign, vcve, lsinw, lcosw in the prior file and use the median e, omegadeg (from the previous .tex file) instead.', ss.logname
      stop
   endif

   if ~finite(ss.planet[i].omega.value) then begin
      printandlog, 'The starting value for omega_' + strtrim(i,2) + ' is not physical. If you are seeing this with an automatically generated prior file after updating EXOFASTv2, delete the the starting points for sign, vcve, lsinw, lcosw in the prior file and use the median e, omegadeg (from the previous .tex file) instead.', ss.logname
      stop
   endif

   ;; translate from e/omega to various e/omega parameterizations
   if ~ss.planet[i].qecosw.userchanged then ss.planet[i].qecosw.value = (ss.planet[i].e.value)^(0.25d0)*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].qesinw.userchanged then ss.planet[i].qesinw.value = (ss.planet[i].e.value)^(0.25d0)*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].secosw.userchanged then ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].sesinw.userchanged then ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].ecosw.userchanged then ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   if ~ss.planet[i].esinw.userchanged then ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   if ~ss.planet[i].vcve.userchanged then ss.planet[i].vcve.value = sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; if the user sets e or omega, 
   ;; we need to make sure our combination of vcve/lcosw/lsinw/lsinw2/sign will recover it
   e1 = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=0B)
   e2 = vcve2e(ss.planet[i].vcve.value,omega=ss.planet[i].omega.value, sign=1B)
  
   if abs(ss.planet[i].e.value - e1) lt 1d-8 then ss.planet[i].sign.value = 0.5 $
   else if abs(ss.planet[i].e.value - e2) lt 1d-8 then ss.planet[i].sign.value = 1.5 $
   else begin
      printandlog, 'ERROR translating starting e and omega to parameterization', ss.logname
      stop
   endelse

   if ~ss.planet[i].lsinw.userchanged and ~ss.planet[i].lcosw.userchanged then begin
      if floor(ss.planet[i].sign.value) then L = sqrt(0.75d0) $
      else L = sqrt(0.25d0)
      ss.planet[i].lsinw.value = L*sin(ss.planet[i].omega.value)
      ss.planet[i].lcosw.value = L*cos(ss.planet[i].omega.value)
   endif

if 0 then begin
   sign = floor(ss.planet[i].sign2.value)
   if ss.planet[i].vcve.value lt 1d0 and ss.planet[i].lsinw.value le 0d0 and sign then begin
      ;; then I need lsinw to be negative and lsinw2 to be positive
      ss.planet[i].lsinw2.value = -ss.planet[i].lsinw.value
      ;; and I need sign to be > 1
      ss.planet[i].sign2.value = 0.5d0
   endif else if ss.planet[i].vcve.value ge 1d0 and ss.planet[i].lsinw.value gt 0d0 and ~sign then begin
      ;; then I need lsinw to be positive and lsinw2 to be negative
      ss.planet[i].lsinw2.value = -ss.planet[i].lsinw.value
      ;; and I need sign to be > 1
      ss.planet[i].sign2.value = 0.5d0
   endif else begin
      ss.planet[i].sign2.value = ss.planet[i].sign.value
      ss.planet[i].lsinw2.value = ss.planet[i].lsinw.value
   endelse
endif else begin
      ss.planet[i].sign2.value = ss.planet[i].sign.value
      ss.planet[i].lsinw2.value = ss.planet[i].lsinw.value
endelse

   ;; how did the user seed rp/rstar?
   if ss.planet[i].p.userchanged then begin
      if ~ss.planet[i].rpearth.userchanged then $
         ss.planet[i].rpearth.value = ss.planet[i].p.value*ss.star[ss.planet[i].starndx].rstar.value/rearth
      ss.planet[i].rpearth.userchanged = 1B
   endif else if ss.planet[i].rpearth.userchanged then begin
      ss.planet[i].p.value = ss.planet[i].rpearth.value*rearth/ss.star[ss.planet[i].starndx].rstar.value
   endif else if ss.planet[i].rp.userchanged then begin
      ss.planet[i].p.value = ss.planet[i].rp.value*rjup/ss.star[ss.planet[i].starndx].rstar.value
      if ~ss.planet[i].rpearth.userchanged then $
         ss.planet[i].rpearth.value = ss.planet[i].p.value*ss.star[ss.planet[i].starndx].rstar.value/rearth
      ss.planet[i].rpearth.userchanged = 1B
   endif else if ss.planet[i].mpearth.userchanged then begin 
      ss.planet[i].rpearth.value = massradius_chen(ss.planet[i].mpearth.value)
      ss.planet[i].p.value = ss.planet[i].rpearth.value*rearth/ss.star[ss.planet[i].starndx].rstar.value
   endif

   ;; how did the user seed the inclination?
   if ss.planet[i].i.userchanged then begin
      sini = sin(ss.planet[i].i.value)
      ss.planet[i].cosi.value = cos(ss.planet[i].i.value)
   endif else if ss.planet[i].ideg.userchanged then begin      
      ss.planet[i].i.value = ss.planet[i].ideg.value*!dpi/180d0     
      sini = sin(ss.planet[i].i.value)
      ss.planet[i].cosi.value = cos(ss.planet[i].i.value)
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
                                       ss.star[ss.planet[i].starndx].mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].k.userchanged then begin
      if ss.planet[i].k.value lt 0d0 then begin
         ss.planet[i].mpsun.value = -ktom2(-ss.planet[i].K.value, ss.planet[i].e.value,$
                                           ss.planet[i].i.value, ss.planet[i].period.value, $
                                           ss.star[ss.planet[i].starndx].mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      endif else begin
         ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                          ss.planet[i].i.value, ss.planet[i].period.value, $
                                          ss.star[ss.planet[i].starndx].mstar.value, GMsun=ss.constants.GMsun/1d6) ;; m_sun
      endelse
      ss.planet[i].logmp.value = alog10(ss.planet[i].mpsun.value)
   endif else if ss.planet[i].rpearth.userchanged then begin
      ss.planet[i].mpearth.value = massradius_chenreverse(ss.planet[i].rpearth.value)
      ss.planet[i].mpsun.value = ss.planet[i].mpearth.value*mearth ;; m_sun
   endif

   ;; how did the user seed a/rstar?
   ss.planet[i].arsun.value=(G*(ss.star[ss.planet[i].starndx].mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                       ;; semi-major axis in r_sun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star[ss.planet[i].starndx].rstar.value    ;; a/rstar (unitless)
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
                                       ss.planet[i].period.value,/tt2tc)

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

;; this is done in exofast_chi2v2.pro
;   if ss.planet[i].fittran and ~ss.noprimary then begin
;      if ~finite(ss.planet[i].chord.value) or $
;         (ss.planet[i].b.value gt (1d0+ss.planet[i].p.value)) or $
;         ~finite(ss.planet[i].cosi.value) then begin
;         printandlog, 'ERROR: the starting values for planet ' + strtrim(i,2) + ' does not transit, but a transit is fit. Revise i, ideg, cosi, b, or chord in the prior file', ss.logname
;         return, 0
;      endif
;   endif

endfor

return, 1

end
