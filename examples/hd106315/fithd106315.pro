pro fithd106315, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=100 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hd106315'])

exofastv2, nplanets=2, tranpath=path + 'n20090101.Kepler.K2.dat',fluxfile=path+'hd106315.flux',$
           priorfile=path + 'hd106315.priors',debug=debug, verbose=verbose, prefix=path + 'hd106315.',$
           maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1],fitrv=[0,0],circular=[0,0],/longcadence

;stop

end
