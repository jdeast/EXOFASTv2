pro fitwasp76, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

;if n_elements(nthin) eq 0 then nthin=50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','wasp76'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','wasp76','fitresults'])

exofastv2, nplanets=1, rvpath=path+'WASP-76b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'wasp76.priors',$
           prefix=outpath+'WASP-76b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,nthread=nthread,$
           /nomist, /torres

end
