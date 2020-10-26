pro fitwasp76, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;if n_elements(nthin) eq 0 then nthin=50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','wasp76'])

exofastv2, nplanets=1, rvpath=path+'WASP-76b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'wasp76.priors',$
           prefix=path+'fitresults' + path_sep()+'WASP-76b.test.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres

end
