pro fitwasp76, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

;if n_elements(nthin) eq 0 then nthin=50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','wasp76'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','wasp76','fitresults'])

;; parallel tempering (ntemps=8) is required because it's bimodal in
;; mass/age, and it can't sample both modes properly without it.
exofastv2, nplanets=1, rvpath=path+'WASP-76b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'wasp76.priors',$
           prefix=outpath+'WASP-76b.PT.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,nthread=nthread,$
           mistsedfile=path+'wasp76.sed',fitthermal=['TESS'],ntemp=8

end
