pro fitkelt6, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=100 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt6'])

exofastv2, nplanets=2, rvpath=path+'KELT-6b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt6.priors',$
           prefix=path+'fitresults' + path_sep()+'KELT-6.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1,1],fittran=[1,0],circular=[0,0], $
           debug=debug,verbose=verbose,$
           fluxfile='kelt6.sed'

end
