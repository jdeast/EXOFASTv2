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
; OPTIONAL INPUTS:
;   STARNDX - The index of the star whose mass is bimodal. Default is
;             0 (the primary). Only relevant for multi-star fits,
;             where MCMCSS.STAR is an NSTARS-element array.
;
; MODIFICATION HISTORY
; 
;  2019/10/11 -- added
;  2026/09/17 -- added multistar support
;-
pro splitpdf, idlfile, masscut, starndx=starndx

restore, idlfile

if n_elements(starndx) eq 0 then starndx = 0

;; NB: mcmcss.star is an NSTARS-element array, so mcmcss.star.mstar.value
;; is [NSTEPS,NSTARS] -- select one star or the cut is meaningless
mstar = mcmcss.star[starndx].mstar.value

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
