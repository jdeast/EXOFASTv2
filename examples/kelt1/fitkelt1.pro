pro fitkelt1, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=10000 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt1'])

exofastv2, nplanets=1, rvpath=path+'KELT-1b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt1.priors',$
           prefix=path+'fitresults' + path_sep()+'KELT-1b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres,fluxfile='kelt1.sed'

end
