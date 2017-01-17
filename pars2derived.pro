pro derivepars, ss

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
AU = 215.094177d0 ;; R_sun

ss.star.mstar.value = 10^ss.star.logmstar.value

for i=0, ss.nplanets-1 do begin
   ss.planet[i].e.value = (ss.planet[i].qecosw.value^2 + ss.planet[i].qesinw.value^2)^2
   ss.planet[i].k.value = 10^ss.planet[i].logk.value
   ss.planet[i].i.value = acos(ss.planet[i].cosi.value)
   ss.planet[i].ideg.value = ss.planet[i].i.value*180d0/!dpi

   ss.planet[i].period.value = 10^ss.planet[i].logp.value
   
   if ss.planet[i].K.value le 0d0 then ss.planet[i].mpsun.value = 0d0 $
   else ss.planet[i].mpsun.value = ktom2(ss.planet[i].K.value, ss.planet[i].e.value,$
                                         ss.planet[i].i.value, ss.planet[i].period.value, $
                                         ss.star.mstar.value)
   
   ss.planet[i].arsun.value=(G*(ss.star.mstar.value+ss.planet[i].mp.value)*ss.planet[i].period.value^2/$
                       (4d0*!dpi^2))^(1d0/3d0)         ;; rsun
   ss.planet[i].ar.value = ss.planet[i].arsun.value/ss.star.rstar.value ;; unitless
   ss.planet[i].a.value = ss.planet[i].arsun.value/AU ;; AU

   
endfor

end
