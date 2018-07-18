pro fitgj9827, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=500 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','gj9827'])

exofastv2, nplanets=3, tranpath=path+'n20090101.Kepler.K2.dat',$
           fluxfile=path+'gj9827.flux',$
           priorfile=path+'gj9827.priors',debug=debug, verbose=verbose, $
           prefix=path+'fitresults' + path_sep()+'GJ9827.',maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1,1],fitrv=[0,0,0],circular=[1,1,1],$
           /longcadence,/noclaret,/earth,/nomist

end
