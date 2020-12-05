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
;mstar = (reform(mcmcss.star.mstar.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]

mstar = mcmcss.star.mstar.value

highmass = where(mstar gt masscut, complement=lowmass)

if highmass[0] eq -1 or lowmass[0] eq -1 then begin
   print, 'masscut (' + strtrim(masscut,2) + ') does not split PDF; median mass is ' + strtrim(median(mstar),2)
   stop
endif

print, 'The probability of the low-mass solution is ' +  strtrim(double(n_elements(lowmass))/n_elements(mstar),2)
print, 'The probability of the high-mass solution is ' +  strtrim(double(n_elements(highmass))/n_elements(mstar),2)

;; high mass solution (cut out low mass solutions)
basename = file_dirname(idlfile) + path_sep() + file_basename(idlfile,'.mcmc.idl') + '.highmass.'
parfile = basename + 'pdf.ps'
covarfile = basename + 'covar.ps'
logname = basename + 'log'
texfile = basename + 'tex'
csvfile = basename + 'csv'
exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile, mask=lowmass
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

;; low mass solution (cut out low mass solutions)
basename = file_dirname(idlfile) + path_sep() + file_basename(idlfile,'.mcmc.idl') + '.lowmass.'
parfile = basename + 'pdf.ps'
covarfile = basename + 'covar.ps'
logname = basename + 'log'
texfile = basename + 'tex'
csvfile = basename + 'csv'
exofast_plotdist_corner, mcmcss, pdfname=parfile, covarname=covarfile,nocovar=nocovar,logname=logname, csvfile=csvfile, mask=highmass
exofast_latextab2, mcmcss, caption=caption, label=label,texfile=texfile

end
