;; flatten a lightcurve using Andrew Vanderberg's Keplerspline
pro flatten, filename,splinespace=splinespace

if n_elements(splinespace) eq 0 then splinespace = 0.75

readcol, filename, time, flux, fluxerr, format='d,d,d'

norm = keplerspline(time, flux, breakp=breakp, ndays=splinespace)

flattened_filename = file_path(filename) + path_sep() + file_basename(filename,'.dat') + '.flattened.dat'

exofast_forprint, time, flux/norm, fluxerr/norm, format='(f0.9,x,f0.9,x,f0.9)',/nocomment, textout=flattened_filename

end
