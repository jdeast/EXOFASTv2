pro fitkelt16, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=10000 ;; 50000
if n_elements(nthin) eq 0 then nthin=4 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt16'])

exofastv2, nplanets=1, rvpath=path+'KELT-16b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt16.priors',$
           prefix=path+'fitresults' + path_sep()+'KELT-16b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres

end
