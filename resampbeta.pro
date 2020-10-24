pro resamp, filename=filename

if n_elements(filename) eq 0 then filename = 'TICID.mcmc.idl'

nsteps=1d3
e = (dindgen(nsteps)+1)/nsteps
b = 3.03d0
a = 0.867d0
beta = e^(a-1d0)*(1d0-e)^(b-1d0)
beta /= beta[0]

plot, e, beta

restore, filename

e = mcmcss.planet[0].e.value
ar = mcmcss.planet[0].ar.value

mine = min(e)
normbeta = mine^(a-1d0)*(1d0-mine)^(b-1d0)
beta = e^(a-1d0)*(1d0-e)^(b-1d0)/normbeta

ran = exofast_random(seed, mcmcss.nsteps)

;; resample for Beta distribution (Kipping, 2013)
goodbeta = where(beta gt ran and e le 1d0-3d0/ar, complement=mask)

basename = file_basename(filename,'.mcmc.idl') + '.beta.'
parfile = basename + 'pdf.ps'
covarfile = basename + 'covar.ps'
logname = basename + 'log'
texfile = basename + 'tex'

exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, angular=angular,csvfile=csvfile, mask=mask
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

stop
end
