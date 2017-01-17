function pars2step, ss

;; derive the stellar radius
ss.star.mstar.value = 10^ss.star.logmstar.value
junk = massradius_yy3(ss.star.mstar.value, ss.star.feh.value, ss.star.age.value, ss.star.teff.value,yyrstar=rstar)
ss.star.rstar.value = rstar
G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010

nplanets = n_elements(ss.planet)
;; derive quantities we'll use later
for i=0, nplanets-1 do begin
   ;; derive quantities we'll use later
   ss.planet[i].logp.value = alog10(ss.planet[i].period.value)
   ss.planet[i].qecosw.value = (ss.planet[i].e.value)^(0.25d0)*cos(ss.planet[i].omega.value)
   ss.planet[i].qesinw.value = (ss.planet[i].e.value)^(0.25d0)*sin(ss.planet[i].omega.value)
   ss.planet[i].logk.value = alog10(ss.planet[i].k.value)

   ;; scale of cosi (scales with projected disk size)
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                    ss.planet[i].i.value, ss.planet[i].period.value, $
                                    ss.star.mstar.value)
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mpsun.value)*ss.planet[i].period.value^2/$
                             (4d0*!dpi^2))^(1d0/3d0)                    ;; rsun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value ;; unitless
   ss.planet[i].a.value = ss.planet[i].arsun.value/215.094177d0         ;; AU
   ss.planet[i].cosi.scale = 1d0/ss.planet[i].ar.value

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
