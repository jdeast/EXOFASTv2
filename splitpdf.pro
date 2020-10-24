;+
; NAME:
;   SPLITPDF
;
; PURPOSE: 
;   Splits an EXOFASTv2 output IDL file and reports two
;   solutions. Intended to summarize bimodal distributions.
;
; CALLING SEQUENCE:
;   splitpdf, idlfile, 1.0
;
; INPUTS:
;   IDLFILE - The name of the EXOFAST output IDL file
;   MASSCUT - The value to cut on, in solar masses
;
; MODIFICATION HISTORY
; 
;  2019/10/11 -- added
;-
pro splitpdf, idlfile, masscut

restore, idlfile

;chi2 = reform(*(mcmcss.chi2),mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
;burnndx = getburnndx(chi2,goodchains=goodchains)
;e = (reform(mcmcss.planet[0].e.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]
;ar = (reform(mcmcss.planet[0].ar.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]

mstar = mcmcss.star.mstar.value
mask1 = where(mstar gt masscut, complement=mask2)

print, 'The probability the low-mass solution is right is ' +  strtrim(double(n_elements(mask2))/n_elements(mstar),2)
print, 'The probability the high-mass solution is right is ' +  strtrim(double(n_elements(mask1))/n_elements(mstar),2)

;; high mass solution
basename = file_basename(idlfile,'.mcmc.idl') + '.highmass.'
parfile = basename + 'pdf.ps'
covarfile = basename + 'covar.ps'
logname = basename + 'log'
texfile = basename + 'tex'
exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile, mask=mask1
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

;; low mass solution
basename = file_basename(idlfile,'.mcmc.idl') + '.lowmass.'
parfile = basename + 'pdf.ps'
covarfile = basename + 'covar.ps'
logname = basename + 'log'
texfile = basename + 'tex'
exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile, mask=mask2
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

end
