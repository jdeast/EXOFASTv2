pro derivepars, ss

;; constants
G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
AU = 215.094177d0 ;; R_sun
mjup = 0.000954638698d0 ;; m_sun
rjup = 0.102792236d0    ;; r_sun
mearth = 0.00000300245 ;; m_sun
rearth = 0.0091705248 ;; r_sun

sigmab = 5.670373d-5/3.839d33*6.9566d10^2 ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
sigmasb = 5.670373d-5

ss.star.mstar.value = 10^ss.star.logmstar.value
ss.star.rhostar.value = ss.star.mstar.value/(ss.star.rstar.value^3)*1.41135837d0
ss.star.logg.value = alog10(G*ss.star.mstar.value/(ss.star.rstar.value^2)*9.31686171d0)
ss.star.lstar.value = 4d0*!dpi*ss.star.rstar.value^2*ss.star.teff.value^4*sigmaB 

;; derive the absolute magnitude and distance
nsteps = n_elements(ss.star.teff.value)
ss.star.parallax.value = 1d3/ss.star.distance.value ;; mas

for i=0, ss.nplanets-1 do begin

   ;; eccentricity and argument of periastron
   if ss.planet[i].qecosw.fit then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
      ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)
   endif else begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endelse
   zero = where(ss.planet[i].e.value eq 0d0,complement=nonzero)
   if zero[0] ne -1 then ss.planet[i].omega.value[zero] = !dpi/2d0 

   ss.planet[i].omegadeg.value = ss.planet[i].omega.value*180d0/!dpi
   ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   ss.planet[i].ecosw.value = ss.planet[i].e.value*cos(ss.planet[i].omega.value)
   ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)
   ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   ss.planet[i].qesinw.value = (ss.planet[i].e.value)^0.25d0*sin(ss.planet[i].omega.value)
   ss.planet[i].qecosw.value = (ss.planet[i].e.value)^0.25d0*cos(ss.planet[i].omega.value)

   ss.planet[i].k.value = 10^ss.planet[i].logk.value
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi
   sini = sin(ss.planet[i].i.value)

   ss.planet[i].period.value = 10^ss.planet[i].logp.value
   
   zero = where(ss.planet[i].logk.value le -6d0,complement=positive)
   if zero[0] ne -1 then begin
      ss.planet[i].mpsun.value[zero] = 0d0
      ss.planet[i].k.value[zero] = 0d0
   endif

   if positive[0] ne -1 then $
      ss.planet[i].mpsun.value[positive] =  ktom2(ss.planet[i].K.value[positive],$
                                                  ss.planet[i].e.value[positive],$
                                                  ss.planet[i].i.value[positive],$
                                                  ss.planet[i].period.value[positive],$
                                                  ss.star.mstar.value[positive])

   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup
   ss.planet[i].msini.value = ss.planet[i].mp.value*sini
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth
   ss.planet[i].msiniearth.value = ss.planet[i].mpearth.value*sini
   ss.planet[i].q.value = ss.planet[i].mpsun.value/ss.star.mstar.value
   
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mp.value)*ss.planet[i].period.value^2/(4d0*!dpi^2))^(1d0/3d0)
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value ;; unitless
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; AU

   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star.rstar.value
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth

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
   ss.planet[i].tp.value = recenter(ss.planet[i].tp.value, medper)
   ss.planet[i].ts.value = recenter(ss.planet[i].ts.value, medper)
   ss.planet[i].ta.value = recenter(ss.planet[i].ta.value, medper)
   ss.planet[i].td.value = recenter(ss.planet[i].td.value, medper)

   ss.planet[i].teq.value = ss.star.teff.value*sqrt(1d0/(2d0*ss.planet[i].ar.value)) ;(f*(1d0-Ab))^(0.25d0)
   ss.planet[i].dr.value = ss.planet[i].ar.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; d/rstar = star-planet separation at transit

   ss.planet[i].fave.value = sigmasb*ss.star.teff.value^4/(ss.planet[i].ar.value*(1d0+ss.planet[i].e.value^2/2d0))^2/1d9    ;; 10^9 erg/s/cm^2

   ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)  ;; eq 7, Winn 2010
   ss.planet[i].bs.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0-ss.planet[i].esinw.value)  ;; eq 8, Winn 2010

   ;; approximate durations taken from Winn 2010 (close enough; these should only be used to schedule observations anyway)
   ss.planet[i].t14.value = ss.planet[i].period.value/!dpi*asin(sqrt((1d0+ss.planet[i].p.value)^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)
   ;; eqs 14, 16, Winn 2010
   t23 = ss.planet[i].period.value/!dpi*asin(sqrt((1d0-ss.planet[i].p.value)^2 - ss.planet[i].b.value^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value)

   ;; no transit, transit duration equation is undefined -- set to zero
   notransit = where(ss.planet[i].b.value gt 1d0+ss.planet[i].p.value)
   if notransit[0] ne -1 then ss.planet[i].t14.value[notransit] = 0d0

   ;; grazing transit, the flat part of transit is zero
   grazing = where(ss.planet[i].b.value gt 1d0-ss.planet[i].p.value)
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
   notransit = where(ss.planet[i].bs.value gt 1d0+ss.planet[i].p.value)
   if notransit[0] ne -1 then ss.planet[i].t14s.value[notransit] = 0d0

   ;; ****************** this is a hack!! *************************
   ;; The above statement should get rid of NaNs, but it doesn't!!
;   bad = where(~finite(ss.planet[i].t14s.value))
;   if bad[0] ne -1 then ss.planet[i].t14s.value[bad] = 0d0
   ;; *************************************************************

   ;; grazing transit, the flat part of transit is zero
   grazing = where(ss.planet[i].bs.value gt 1d0-ss.planet[i].p.value)
   if grazing[0] ne -1 then t23s[grazing] = 0d0

   ss.planet[i].taus.value = (ss.planet[i].t14s.value-t23s)/2d0
   ss.planet[i].tfwhms.value = ss.planet[i].t14s.value-ss.planet[i].taus.value

   if ss.planet[i].tfwhms.derive and (where(~finite(ss.planet[i].tfwhms.value)))[0] ne -1 then stop
   if ss.planet[i].taus.derive and (where(~finite(ss.planet[i].taus.value)))[0] ne -1 then stop

   ss.planet[i].psg.value = (ss.star.rstar.value+ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2) ;; eq 9, Winn 2010
   ss.planet[i].ps.value = (ss.star.rstar.value-ss.planet[i].rpsun.value)/ss.planet[i].arsun.value*(1d0 - ss.planet[i].esinw.value)/(1d0-ss.planet[i].e.value^2)  ;; eq 9, Winn 2010

   ss.planet[i].rhop.value = ss.planet[i].mpsun.value/(ss.planet[i].rpsun.value^3)*1.41135837d0
   ss.planet[i].loggp.value = alog10(G*ss.planet[i].mp.value/ss.planet[i].rp.value^2*9.31686171d0) ;; cgs
   ss.planet[i].safronov.value = ss.planet[i].ar.value*ss.planet[i].q.value/ss.planet[i].p.value

   ;; depth != delta if grazing (ignore limb darkening)
   ss.planet[i].delta.value = ss.planet[i].p.value^2
   for j=0, ss.nsteps-1 do begin
      exofast_occultquad, ss.planet[i].b.value[j], 0d0, 0d0, ss.planet[i].p.value[j],mu1
      ss.planet[i].depth.value[j] = 1d0-mu1[0]
   endfor

   if ss.planet[i].fittran then begin
      

   endif


endfor

end
