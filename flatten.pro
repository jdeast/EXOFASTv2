pro flatten, filename

readcol, filename, time, flux, format='d,d',delimiter=',',skip=1

flatflux = flux/median(flux,48)

forprint, time+2454833d0, flatflux, flatflux*0d0 + stddev(flatflux), format='(f0.9,x,f0.9,x,f0.9)',/nocomment, textout=filename + '.dat3'

end
