pro rhostar2ar, rhostar, period, urhostar=urhostar, p=p, rhop=rhop

;; rhostar is in g/cm^3

constants = mkconstants()

if n_elements(p) eq 0 then p = rhostar*0d0
if n_elements(rhop) eq 0 then rhop = rhostar*0d0


rhostararr = rhostar + randomn(seed,1d6)*urhostar + p^3*rhop
ar = (rhostararr*(constants.G*(period*86400d0)^2)/(3d0*!dpi))^(1d0/3d0)

npoints = 1d6
ar = 11.24d0 + randomn(seed,npoints)*0.05d0
period = 1.46712158 + randomn(seed,npoints)*0.00000049d0
rhostar2 = ar^3*3d0*!dpi/(constants.G*(period*86400d0)^2) ;; cgs


rstar = 0.328d0 + randomn(seed,npoints)*0.011d0
mstar = 0.323d0 + randomn(seed,npoints)*0.015d0
rhostar3 = mstar/rstar^3*constants.rhosun
ar2 = (rhostar3*(constants.G*(period*86400d0)^2)/(3d0*!dpi))^(1d0/3d0)

print, median(ar), stdev(ar)
print, median(ar2), stdev(ar2)

print, median(rhostar2), stdev(rhostar2)
print, median(rhostar3), stdev(rhostar3)


stop

end
