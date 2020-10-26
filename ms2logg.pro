function ms2logg, mstar, rstar

constants = mkconstants()
return, alog10(mstar/rstar^2*constants.gravitysun)

end
