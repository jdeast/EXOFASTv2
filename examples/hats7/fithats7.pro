pro fithats7, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=100 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hats7'])

exofastv2, nplanets=1, rvpath=path+'HATS-7b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'hats7.priors',$
           prefix=path+'fitresults' + path_sep()+'HATS-7b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose;,$
;           /nomist, /torres

end
