pro fitkelt16, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt16'])

if n_elements(outpath) eq 0 then $
   outpath = 'modeling' + path_sep() + 'kelt16' + path_sep() + 'fitresults' + path_sep()

;; generally i wouldn't recommend using /nomist and /torres
exofastv2, nplanets=1, rvpath=path+'KELT-16b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt16.priors',$
           prefix=outpath+'KELT-16b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,nthread=nthread,$
           /nomist, /torres

end
