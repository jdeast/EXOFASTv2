pro fitkelt16, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt16'])

;; generally i wouldn't recommend using /nomist and /torres
exofastv2, nplanets=1, rvpath=path+'KELT-16b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt16.priors',$
           prefix=path+'fitresults' + path_sep()+'KELT-16b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,nthread=nthread,$
           /nomist, /torres

end
