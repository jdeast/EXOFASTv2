pro fitkelt12, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, mist=mist, gaia=gaia

;; default to a very short run (not well-mixed or publication quality)                         
if n_elements(maxsteps) eq 0 then maxsteps=2500 ;; 50000                                      
if n_elements(nthin) eq 0 then nthin=20 ;; 50                                                  
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt12'])

exofastv2, nplanets=1, rvpath=path + 'KELT-12b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt12.priors',$
           prefix=path + 'fitresults/KELT-12.',$
           maxsteps=maxsteps,nthin=nthin,circular=[1],$
           debug=debug, verbose=verbose,$
           mistsedfile=path+'kelt12.sed',/fitslope,/ttv

end
