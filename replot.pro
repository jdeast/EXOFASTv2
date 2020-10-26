pro replot, savefile

savefile = 'HAT-P-3b.torres.mcmc.idl'
savefile = 'HAT-3b.YY.mcmc.idl'
;savefile = 'HAT-P-3b.MIST.mcmc.idl'
restore, savefile

;prefix = mcmcss.prefix
prefix = 'test.'
basename = file_basename(prefix)

;; output filenames
label = "tab:" + basename
caption = "Median values and 68\% confidence interval for " + basename
parfile = prefix + 'pdf.ps'
covarfile = prefix + 'covar.ps'
chainfile = prefix + 'chain.ps'
texfile = prefix + 'median.tex'

exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, angular=angular

spawn, 'gv ' + covarfile + ' &'

stop
end
