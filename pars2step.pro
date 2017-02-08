function pars2step, ss

;; derive the stellar mass
ss.star.mstar.value = 10^ss.star.logmstar.value
if ss.star.mstar.scale ne 0d0 then $
   ss.star.logmstar.scale = ss.star.mstar.scale/(alog(10d0)*ss.star.mstar.value)

;; if a starting value for the stellar radius was given, 
;; use it to derive a starting point for the age
if ss.star.rstar.value ne 0d0 then begin
   ntries = 100
   age = dindgen(ntries)/(ntries-1)*13.82
   rstars = dblarr(ntries)
   for i=0, ntries-1 do begin
      junk = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, age[i], ss.star.teff.value,yyrstar=rstar)
      rstars[i] = rstar
   endfor
   junk = min(abs(rstars-ss.star.rstar.value),ndx)
   ss.star.age.value = age[ndx]
endif else begin
   ;; otherwise derive the stellar radius
   junk = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, ss.star.age.value, ss.star.teff.value,yyrstar=rstar)
   ss.star.rstar.value = rstar
endelse

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010

nplanets = n_elements(ss.planet)
;; derive quantities we'll use later
for i=0, nplanets-1 do begin

   ;; derive quantities we'll use later
   if ss.planet[i].logp.value eq 0d0 and ss.planet[i].period.value gt 0d0 then begin
      ss.planet[i].logp.value = alog10(ss.planet[i].period.value)

      ;; translate the stepping scale from P to logP
      if ss.planet[i].period.scale ne 0d0 then $
         ss.planet[i].logp.scale = ss.planet[i].period.scale/(alog(10d0)*ss.planet[i].period.value)

   endif

   if ss.planet[i].logk.value eq 1d0 and ss.planet[i].k.value gt 0d0 then begin
      ss.planet[i].logk.value = alog10(ss.planet[i].k.value)

      ;; translate the stepping scale from K to logK
      if ss.planet[i].k.scale ne 0d0 then $
         ss.planet[i].logk.scale = ss.planet[i].k.scale/(alog(10d0)*ss.planet[i].k.value)
   endif

   ;; translate from omegadeg to omega
   if ss.planet[i].omegadeg.value ne 0 then begin
      ss.planet[i].omega.value = ss.planet[i].omegadeg.value*!pi/180d0
      if ss.planet[i].omegadeg.scale ne 0 then ss.planet[i].omega.scale = ss.planet[i].omegadeg.scale*!pi/180d0
   endif

   if ss.planet[i].qecosw.value eq 0d0 then ss.planet[i].qecosw.value = (ss.planet[i].e.value)^(0.25d0)*cos(ss.planet[i].omega.value)
   if ss.planet[i].qesinw.value eq 0d0 then ss.planet[i].qesinw.value = (ss.planet[i].e.value)^(0.25d0)*sin(ss.planet[i].omega.value)
   if ss.planet[i].secosw.value eq 0d0 then ss.planet[i].secosw.value = sqrt(ss.planet[i].e.value)*cos(ss.planet[i].omega.value)
   if ss.planet[i].sesinw.value eq 0d0 then ss.planet[i].sesinw.value = sqrt(ss.planet[i].e.value)*sin(ss.planet[i].omega.value)

   ;; scale of cosi (scales with projected disk size)
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                                                           ss.planet[i].i.value, ss.planet[i].period.value, $
                                                                           ss.star.mstar.value)
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                                                                    (4d0*!dpi^2))^(1d0/3d0)                    ;; rsun
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value                                        ;; unitless
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].a.value = ss.planet[i].arsun.value/215.094177d0         ;; AU
   if ss.planet[i].logp.value eq 0d0 then ss.planet[i].cosi.scale = 1d0/ss.planet[i].ar.value

   if ss.planet[i].tc.value eq 0d0 and ss.planet[i].tp.value ne 0 then begin
      phase = exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value, /primary)
      ss.planet[i].tc.value = ss.planet[i].tp.value + ss.planet[i].period.value*phase
   endif

   ;; limit eccentricity to avoid collision with star during periastron
   ;; the ignored tidal effects would become important long before this,
   ;; but this prevents numerical problems compared to the e < 1 constraint
   ;; the not/lt (instead of ge) robustly handles NaNs too
   ;; abs(p) because p is allowed to be negative to eliminate bias
   ;; if not (ss.planet[i].e.value lt (1d0-1d0/ss.planet[i].ar.value-abs(ss.planet[i].p.value)/$
   ;;                            ss.planet[i].ar.value)) then return, -1

   ;; ;; tighter (empirical) constraint on eccentricity (see Eastman 2013)
   ;; if ss.tides and ss.planet[i].e.value gt (1d0-3d0/ss.planet[i].ar.value) then begin
   ;;    ss.planet[i].e.value = 0d0
   ;;    ss.planet[i].qesinw.value = 0d0
   ;;    ss.planet[i].qecosw.value = 0d0
   ;; endif

   ;; if ss.planet[i].e.value eq 0d0 then ss.planet[i].omega.value = !dpi/2d0 $
   ;; else ss.planet[i].omega.value = atan(ss.planet[i].qesinw.value,ss.planet[i].qecosw.value)

   ;; ;; time of periastron
   ;; ss.planet[i].phase.value=exofast_getphase(ss.planet[i].e.value,ss.planet[i].omega.value,/pri)
   ;; ss.planet[i].tp.value = ss.planet[i].tc.value - ss.planet[i].period.value*ss.planet[i].phase.value

endfor

end
