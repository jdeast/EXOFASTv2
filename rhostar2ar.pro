pro rhostar2ar, rhostar, period, urhostar=urhostar, p=p, rhop=rhop

;; rhostar is in g/cm^3

constants = mkconstants()

rhostararr = rhostar + randomn(seed,1d6)*urhostar + p^3*rhop
ar = (rhostararr*(constants.G*(period*86400d0)^2)/(3d0*!dpi))^(1d0/3d0)

stop

end
