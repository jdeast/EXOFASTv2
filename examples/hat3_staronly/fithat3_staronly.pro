pro fithat3_staronly, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
             nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3_staronly'])

if n_elements(outpath) eq 0 then $
   outpath = 'modeling' + path_sep() + 'hat3_staronly' + path_sep() + 'fitresults' + path_sep()

;; the SED and prior files were created with this command:
;; we won't do it here so as not to pollute the EXOFASTv2 repo
;mkticsed, '311035838', priorfile=path+'HAT-3.priors', sedfile=path+'HAT-3.sed'

exofastv2, nplanets=0, priorfile=path+'HAT-3.priors', mistsedfile=path+'HAT-3.sed', $
           prefix=outpath +'HAT-3.staronly.', maxsteps=maxsteps, ntemp=8, $
           nthin=nthin, debug=debug, verbose=verbose, nthread=nthread

end
