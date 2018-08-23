pro derivepars, ss

ss.star.mstar.value = 10^ss.star.logmstar.value
ss.star.rhostar.value = ss.star.mstar.value/(ss.star.rstar.value^3)*ss.constants.rhosun ;; rho_sun
ss.star.logg.value = alog10(ss.star.mstar.value/(ss.star.rstar.value^2)*ss.constants.gravitysun) ;; cgs
ss.star.lstar.value = 4d0*!dpi*ss.star.rstar.value^2*ss.star.teff.value^4*ss.constants.sigmab/ss.constants.lsun*ss.constants.rsun^2 ;; lsun
;ss.star.fbol.value = (ss.star.lstar.value/ss.constants.lsun)/(4d0*!dpi*(ss.star.distance.value/ss.constants.pc)^2) ;; cgs

;; derive the age
if ss.mist and not ss.star.age.fit then begin
   for i=0L, ss.nsteps-1 do begin
      mistchi2 = massradius_mist(ss.star.eep.value[i],ss.star.mstar.value[i],$
                                 ss.star.initfeh.value[i],ss.star.age.value[i],$
                                 ss.star.teff.value[i],ss.star.rstar.value[i],$
                                 ss.star.feh.value[i],mistage=mistage,fitage=ss.star.age.fit)
      ss.star.age.value[i] = mistage
   endfor
endif

;; derive the absolute magnitude and distance
if ss.star.distance.fit then ss.star.parallax.value = 1d3/ss.star.distance.value ;; mas
if ss.star.parallax.fit then ss.star.distance.value = 1d3/ss.star.parallax.value ;; pc

for j=0, ss.ntel-1 do begin
   positive = where(ss.telescope[j].jittervar.value gt 0d0)
   ss.telescope[j].jitter.value[positive] = sqrt(ss.telescope[j].jittervar.value[positive])
endfor

;; find the extent of the data
minbjd = !values.d_infinity
maxbjd = -!values.d_infinity
for i=0, ss.ntel-1 do begin
   tmpminbjd = min((*ss.telescope[i].rvptrs).bjd,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor
for i=0, ss.ntran-1 do begin
   tmpminbjd = min((*ss.transit[i].transitptrs).bjd,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor
for i=0, ss.nastrom-1 do begin
   tmpminbjd = min((*ss.astrom[i].astromptrs).bjdtdb,max=tmpmaxbjd)
   if tmpminbjd lt minbjd then minbjd = tmpminbjd
   if tmpmaxbjd gt maxbjd then maxbjd = tmpmaxbjd
endfor

starbb36 = exofast_blackbody(ss.star.teff.value,replicate(3550d0/1d9,ss.nsteps),/wave)
starbb45 = exofast_blackbody(ss.star.teff.value,replicate(4493d0/1d9,ss.nsteps),/wave)

for i=0, ss.nplanets-1 do begin

   ;; eccentricity and argument of periastron
   if ss.planet[i].secosw.fit and ss.planet[i].sesinw.fit then begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endif else if ss.planet[i].omega.fit and ss.planet[i].vvcirc.fit then begin
      a = ss.planet[i].vvcirc.value^2*sin(ss.planet[i].omega.value)^2 + 1d0
      b = 2d0*ss.planet[i].vvcirc.value^2*sin(ss.planet[i].omega.value)
      c = ss.planet[i].vvcirc.value^2-1d0

      e1 = (-b + sqrt(b^2 - 4d0*a*c))/(2d0*a)
      e2 = (-b + sqrt(b^2 - 4d0*a*c))/(2d0*a)

      for j=0L, n_elements(a)-1 do begin
         

         if ~finite(e1[j]) or e1[j] lt 0 then begin
            ss.planet[i].e.value[j] = e2[j] ;; e1 bad, e2 good; use e2
         endif else begin
            if ~finite(e2[j]) or e2[j] lt 0 then begin
               ss.planet[i].e.value[j] = e1[j] ;; e2 bad, e1 good, use e1
            endif else begin
               ;; both e1 and e2 are good
               ;; use a reproduceable seed to get a 50/50 chance 
               ;; to recreate this choice later
               ;; sort of an abuse of RNG seeds, but not bad...
               random = exofast_random((ss.planet[i].vvcirc.value[j]-floor(ss.planet[i].vvcirc.value[j]))*(2ULL^64-1))
               if random ge 0.5 then ss.planet[i].e.value[j] = e1[j] $
               else ss.planet[i].e.value[j] = e2[j]
            endelse
         endelse
      endfor
   endif else if ss.planet[i].qecosw.fit and ss.planet[i].qesinw.fit then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
      ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)
   endif

   zero = where(ss.planet[i].e.value eq 0d0,complement=nonzero)
   if zero[0] ne -1 then ss.planet[i].omega.value[zero] = !dpi/2d0 

   ss.planet[i].lambdadeg.value = ss.planet[i].lambda.value*180d0/!dpi
   ss.planet[i].omegadeg.value = ss.planet[i].omega.value*180d0/!dpi
   ss.planet[i].bigomegadeg.value = ss.planet[i].bigomega.value*180d0/!dpi
   ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)
   ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   ss.planet[i].qesinw.value = (ss.planet[i].e.value)^0.25d0*sin(ss.planet[i].omega.value)
   ss.planet[i].qecosw.value = (ss.planet[i].e.value)^0.25d0*cos(ss.planet[i].omega.value)

   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi
   sini = sin(ss.planet[i].i.value)

   ss.planet[i].period.value = 10^ss.planet[i].logp.value
   
   if ss.planet[i].logk.fit then ss.planet[i].k.value = 10^ss.planet[i].logk.value $
   else ss.planet[i].logk.value = alog10(ss.planet[i].k.value)

   zero = where(ss.planet[i].k.value eq 0d0)
   if zero[0] ne -1 then begin
      ss.planet[i].mpsun.value[zero] = 0d0
   endif

   positive = where(ss.planet[i].k.value gt 0)
   if positive[0] ne -1 then $
      ss.planet[i].mpsun.value[positive] =  ktom2(ss.planet[i].K.value[positive],$
                                                  ss.planet[i].e.value[positive],$
                                                  ss.planet[i].i.value[positive],$
                                                  ss.planet[i].period.value[positive],$
                                                  ss.star.mstar.value[positive], GMsun=ss.constants.GMsun/1d6) ;; convert G to SI

   ;; guard against the Lucy-Sweeney type bias by allowing negative K/masses
   negative = where(ss.planet[i].k.value lt 0)
   if negative[0] ne -1 then $
      ss.planet[i].mpsun.value[negative] =  -ktom2(-ss.planet[i].K.value[negative],$
                                                   ss.planet[i].e.value[negative],$
                                                   ss.planet[i].i.value[negative],$
                                                   ss.planet[i].period.value[negative],$
                                                   ss.star.mstar.value[negative], GMsun=ss.constants.GMsun/1d6) ;; convert G to SI

   ss.planet[i].mp.value = ss.planet[i].mpsun.value/ss.constants.GMJupiter*ss.constants.GMSun    ;; m_jupiter
   ss.planet[i].msini.value = ss.planet[i].mp.value*sini                                         ;; m_jupiter
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/ss.constants.GMearth*ss.constants.GMSun ;; m_earth
   ss.planet[i].msiniearth.value = ss.planet[i].mpearth.value*sini                               ;; m_earth
   ss.planet[i].q.value = ss.planet[i].mpsun.value/ss.star.mstar.value                           ;; unitless
   
   G = ss.constants.GMsun/ss.constants.rsun^3*ss.constants.day^2
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/(4d0*!dpi^2))^(1d0/3d0) ;; (a1 + a2)/rsun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value                                                        ;; (a1 + a2)/rstar
   ss.planet[i].a.value = ss.planet[i].arsun.value*ss.constants.rsun/ss.constants.au                                           ;; AU
   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star.rstar.value                                                         ;; r_sun
   ss.planet[i].rp.value = ss.planet[i].rpsun.value*ss.constants.Rsun/ss.constants.RJupiter                                    ;; r_jupiter
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value*ss.constants.Rsun/ss.constants.REarth                                 ;; r_earth

   ;; time of periastron
   ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)  
   phase2 = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/secondary)
   phasea = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/ascending)
   phased = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/descending)

   ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value)
   ss.planet[i].ts.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phase2)
   ss.planet[i].ta.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phasea)
   ss.planet[i].td.value = ss.planet[i].tc.value - ss.planet[i].period.value*(ss.planet[i].phase.value-phased)

   ;; it's possible tp,ts,ta,td could be split down the middle 
   ;; then the median would be meaningless -- correct that
   medper = median(ss.planet[i].period.value)
   ss.planet[i].tp.value = exofast_recenter(ss.planet[i].tp.value, medper)
   ss.planet[i].ts.value = exofast_recenter(ss.planet[i].ts.value, medper)
   ss.planet[i].ta.value = exofast_recenter(ss.planet[i].ta.value, medper)
   ss.planet[i].td.value = exofast_recenter(ss.planet[i].td.value, medper)

   ss.planet[i].teq.value = ss.star.teff.value*sqrt(1d0/(2d0*ss.planet[i].ar.value)) ;(f*(1d0-Ab))^(0.25d0)
   ss.planet[i].dr.value = ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; d/rstar = star-planet separation at transit

   ss.planet[i].fave.value = ss.constants.sigmab*ss.star.teff.value^4/(ss.planet[i].ar.value*(1d0+ss.planet[i].e.value^2/2d0))^2/1d9    ;; 10^9 erg/s/cm^2

   ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)  ;; eq 7, Winn 2010
   ss.planet[i].bs.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)  ;; eq 8, Winn 2010

   ;; approximate durations taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].t14.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+ss.planet[i].p.value)^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)
   ;; eqs 14, 16, Winn 2010
   t23 = ss.planet[i].period.value/!dpi*asin(sqrt((1d0-ss.planet[i].p.value)^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; no transit, transit duration equation is undefined -- set to zero
   notransit = where(abs(ss.planet[i].b.value) gt 1d0+abs(ss.planet[i].p.value))
   if notransit[0] ne -1 then ss.planet[i].t14.value[notransit] = 0d0

   ;; grazing transit, the flat part of transit is zero
   grazing = where(abs(ss.planet[i].b.value) gt 1d0-abs(ss.planet[i].p.value))
   if grazing[0] ne -1 then t23[grazing] = 0d0

   ss.planet[i].tau.value = (ss.planet[i].t14.value-t23)/2d0
   ss.planet[i].tfwhm.value = ss.planet[i].t14.value-ss.planet[i].tau.value

 ;  ;; hack....
 ;  if (where(~finite(ss.planet[i].tfwhm.value)))[0] ne -1 then stop

   ss.planet[i].ptg.value = (ss.star.rstar.value+ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 + ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2) ;; eq 9, Winn 2010
   ss.planet[i].pt.value = (ss.star.rstar.value-ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 + ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2)  ;; eq 9, Winn 2010

   ;; approximate durations taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].t14s.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+ss.planet[i].p.value)^2 - ss.planet[i].bs.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value) ;; eqs 14, 16, Winn 2010
   t23s = ss.planet[i].period.value/!dpi*asin(sqrt((1d0-ss.planet[i].p.value)^2 - ss.planet[i].bs.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)   

   ;; no transit => transit duration equation is undefined -- set to zero
   notransit = where(abs(ss.planet[i].bs.value) gt 1d0+abs(ss.planet[i].p.value))
   if notransit[0] ne -1 then ss.planet[i].t14s.value[notransit] = 0d0

   ;; ****************** this is a hack!! *************************
   ;; The above statement should get rid of NaNs, but it doesn't!!
;   bad = where(~finite(ss.planet[i].t14s.value))
;   if bad[0] ne -1 then ss.planet[i].t14s.value[bad] = 0d0
   ;; *************************************************************

   ;; grazing transit, the flat part of transit is zero
   grazing = where(abs(ss.planet[i].bs.value) gt 1d0-abs(ss.planet[i].p.value))
   if grazing[0] ne -1 then t23s[grazing] = 0d0

   ss.planet[i].taus.value = (ss.planet[i].t14s.value-t23s)/2d0
   ss.planet[i].tfwhms.value = ss.planet[i].t14s.value-ss.planet[i].taus.value

   if ss.planet[i].tfwhms.derive and (where(~finite(ss.planet[i].tfwhms.value)))[0] ne -1 then stop
   if ss.planet[i].taus.derive and (where(~finite(ss.planet[i].taus.value)))[0] ne -1 then stop

   ss.planet[i].psg.value = (ss.star.rstar.value+ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2) ;; eq 9, Winn 2010
   ss.planet[i].ps.value = (ss.star.rstar.value-ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2)  ;; eq 9, Winn 2010

   ss.planet[i].rhop.value = ss.planet[i].mpsun.value/(ss.planet[i].rpsun.value^3)*ss.constants.rhosun ;; cgs
   ss.planet[i].loggp.value = alog10(ss.planet[i].mpsun.value/ss.planet[i].rpsun.value^2*ss.constants.gravitysun) ;; cgs
   ss.planet[i].safronov.value = ss.planet[i].ar.value*ss.planet[i].q.value/ss.planet[i].p.value ;; unitless

   ;; depth != delta if grazing (ignore limb darkening)
   ss.planet[i].delta.value = ss.planet[i].p.value^2
   for j=0L, ss.nsteps-1L do begin
      exofast_occultquad, abs(ss.planet[i].b.value[j]), 0d0, 0d0, ss.planet[i].p.value[j],mu1
      ss.planet[i].depth.value[j] = 1d0-mu1[0]
   endfor

   ;; find the ideal Tc that minimizes the Tc uncertainty and the
   ;; covariance between Tc and Period
   minepoch = floor((minbjd - median(ss.planet[i].tc.value))/median(ss.planet[i].period.value))-1
   maxepoch =  ceil((maxbjd - median(ss.planet[i].tc.value))/median(ss.planet[i].period.value))+1
   mincovar = !values.d_infinity
   for epoch=minepoch, maxepoch do begin
      corr = correlate(transpose([[ss.planet[i].tc.value+epoch*ss.planet[i].period.value],[ss.planet[i].period.value]]))
      if abs(corr[0,1]) lt mincovar then begin
         mincovar = abs(corr[0,1])
         bestepoch = epoch
      endif
   endfor
   ;; apply the best
   ss.planet[i].t0.value = ss.planet[i].tc.value + bestepoch*ss.planet[i].period.value

   ;; blackbody eclipse depths
   planetbb36 = exofast_blackbody(ss.planet[i].teq.value,replicate(3550d0/1d9,ss.nsteps),/wave)
   planetbb45 = exofast_blackbody(ss.planet[i].teq.value,replicate(4493d0/1d9,ss.nsteps),/wave)
   x = ss.planet[i].p.value^2*planetbb36/starbb36
   ss.planet[i].eclipsedepth36.value = x/(1d0+x)*1d6
   x = ss.planet[i].p.value^2*planetbb45/starbb45
   ss.planet[i].eclipsedepth45.value = x/(1d0+x)*1d6
   
;   if ss.planet[i].fittran then begin
;      
;
;   endif


endfor

end
