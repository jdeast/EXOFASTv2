pro recalc, idlfile

if n_elements(idlfile) eq 0 then $;idlfile = 'EPIC229021605.c.mcmc.idl'
   idlfile = 'GJ1132.c.noteffprior.demc2.mcmc.idl'
restore, idlfile
prefix='EPIC229021605.c2.'
;; output filenames
basename = file_basename(prefix)
label = "tab:" + basename
caption = "Median values and 68\% confidence interval for " + basename
parfile = prefix + 'pdf.ps'
covarfile = prefix + 'covar.ps'
texfile = prefix + 'median.tex'
chainfile = prefix + 'chain.ps'

if keyword_set(covar) then nocovar = 0 $
else nocovar = 1

exofast_plotdist2, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile
exofast_plotchains, mcmcss, chainfile=chainfile

end
