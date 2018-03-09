pro fitkelt4, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=50000 ;; 50000
if n_elements(nthin) eq 0 then nthin=20 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt4rvonly'])

exofastv2, nplanets=1, rvpath=path+'KELT-4b.*.rv',$
           priorfile=path+'kelt4.priors',$
           prefix=path+'KELT-4.MIST.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /noyy, /mist

end
