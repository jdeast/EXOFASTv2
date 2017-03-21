function step2pars, ss, verbose=verbose

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
AU = 215.094177d0 ;; R_sun
mjup = 0.000954638698d0 ;; m_sun
rjup = 0.102792236d0    ;; r_sun
mearth = 0.00000300245 ;; m_sun
rearth = 0.0091705248 ;; r_sun
mindist = dblarr(ss.nplanets)
maxdist = dblarr(ss.nplanets)

;; derive stellar parameters
ss.star.mstar.value = 10^ss.star.logmstar.value
ss.star.logg.value = alog10(G*ss.star.mstar.value/(ss.star.rstar.value^2)*9.31686171d0)
;; derive the distance from lstar
sigmab = 5.670373d-5/3.839d33*6.9566d10^2 ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
ss.star.lstar.value = 4d0*!dpi*ss.star.rstar.value^2*ss.star.teff.value^4*sigmaB    ;; L_sun
ss.star.parallax.value = 1d3/ss.star.distance.value ;; mas

for i=0, ss.nplanets-1 do begin

   ss.planet[i].k.value = 10^ss.planet[i].logk.value
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].period.value = 10^ss.planet[i].logp.value

   ;; bound tc to be within 1 period of the starting value for tc
   if ss.planet[i].tc.prior ne 0d0 then begin
      nper = round((ss.planet[i].tc.prior - ss.planet[i].tc.value)/ss.planet[i].period.value)
      ss.planet[i].tc.value -= nper*ss.planet[i].period.value
   endif

   ;; derive quantities we'll use later
   if ss.planet[i].qecosw.fit then begin
      ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2 
      ss.planet[i].omega.value = atan(ss.planet[i].qecosw.value,ss.planet[i].qesinw.value)
   endif else begin
      ss.planet[i].e.value = ss.planet[i].secosw.value^2 + ss.planet[i].sesinw.value^2
      ss.planet[i].omega.value = atan(ss.planet[i].sesinw.value,ss.planet[i].secosw.value)
   endelse

   if ss.planet[i].k.value le 0d0 then ss.planet[i].mpsun.value = 0d0 $
   else ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                         ss.planet[i].i.value, ss.planet[i].period.value, $
                                         ss.star.mstar.value)

;   if ss.planet[i].mpsun.value gt 0.08d0 then begin
;      if keyword_set(verbose) then print, 'Planet ' + strtrim(i,2) + ' mass (' + strtrim(ss.planet[i].mpsun.value/mjup,2) + ' M_J) above hydrogen burning limit'
;      return, -1 ;; planet above the hydrogen burning limit
;   endif

   ss.planet[i].mp.value = ss.planet[i].mpsun.value/mjup
   ss.planet[i].mpearth.value = ss.planet[i].mpsun.value/mearth

   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mp.value)*ss.planet[i].period.value^2/$
                       (4d0*!dpi^2))^(1d0/3d0)         ;; rsun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value ;; unitless
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; AU

   ss.planet[i].rpsun.value = ss.planet[i].p.value*ss.star.rstar.value
   ss.planet[i].rp.value = ss.planet[i].rpsun.value/rjup
   ss.planet[i].rpearth.value = ss.planet[i].rpsun.value/rearth

   ss.planet[i].esinw.value = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
   ss.planet[i].b.value = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+ss.planet[i].esinw.value) ;; eq 7, Winn 2010

   ;; no transit and we're fitting a transit! This causes major problems; exclude this model
   if ss.planet[i].fittran and (ss.planet[i].b.value gt (1d0+ss.planet[i].p.value)) then begin
      if keyword_set(verbose) then print, 'Planet does not transit!'
      return, -1
   endif

   ;; limit eccentricity to avoid collision with star during periastron
   ;; the ignored tidal effects would become important long before this,
   ;; but this prevents numerical problems compared to the e < 1 constraint
   ;; the not/lt (instead of ge) robustly handles NaNs too
   ;; abs(p) because p is allowed to be negative to eliminate bias
   if not (ss.planet[i].e.value lt (1d0-1d0/ss.planet[i].ar.value-abs(ss.planet[i].p.value)/$
                              ss.planet[i].ar.value)) then begin
      if keyword_set(verbose) then print, 'Planet ' + strtrim(i,2) + ' will collide with the star! e=' + strtrim(ss.planet[i].e.value,2) + '; a/Rstar=' + strtrim(ss.planet[i].ar.value,2) + '; Rp/Rstar=' + strtrim(ss.planet[i].p.value,2) + '; Rstar=' + strtrim(ss.star.rstar.value,2)
      if keyword_set(verbose) then print, 'Rstar is derived from YY isochrones; adjust starting values for Age, Teff, [Fe/H], or Mstar to change it'
      return, -1
   endif

   ;; for multi-planet systems, make sure they don't enter each other's hill spheres
   ;; if mp unknown, mp=0 => hill radius=0 => planets can't cross orbits
   ;; **** ignores mutual inclination; a priori excludes systems like Neptune and Pluto!! ****
   hillradius = (1d0-ss.planet[i].e.value)*ss.planet[i].a.value*(ss.planet[i].mpsun.value/(3d0*ss.star.mstar.value))^(1d0/3d0)
   mindist[i] = (1d0-ss.planet[i].e.value)*ss.planet[i].a.value - hillradius
   maxdist[i] = (1d0+ss.planet[i].e.value)*ss.planet[i].a.value + hillradius
   for j=i-1,0,-1 do begin
      if ((mindist[i] ge mindist[j]) and (mindist[i] le maxdist[j])) or $
         ((maxdist[i] ge mindist[j]) and (maxdist[i] le maxdist[j])) or $
         ((mindist[j] ge mindist[i]) and (mindist[j] le maxdist[i])) or $
         ((maxdist[j] ge mindist[i]) and (maxdist[j] le maxdist[i])) then begin
         if keyword_set(verbose) then print, 'Planets ' + strtrim(j,2) + ' (' + strtrim(mindist[j],2) + ',' + strtrim(maxdist[j],2) + ') and ' + strtrim(i,2) + ' (' + strtrim(mindist[i],2) + ',' + strtrim(maxdist[i],2) + ') cross paths'
         return, -1
      endif
   endfor

   ;; tighter (empirical) constraint on eccentricity (see Eastman 2013)
   if ss.tides and ss.planet[i].e.value gt (1d0-3d0/ss.planet[i].ar.value) then begin
      ss.planet[i].e.value = 0d0
      ss.planet[i].qesinw.value = 0d0
      ss.planet[i].qecosw.value = 0d0
   endif

   if ss.planet[i].e.value eq 0d0 then ss.planet[i].omega.value = !dpi/2d0

   ;; time of periastron
   ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)
   ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*ss.planet[i].phase.value

endfor

return, 1
end
