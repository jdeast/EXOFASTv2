pro fitkelt12, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt12'])

exofastv2, nplanets=1, rvpath=path + 'KELT-12b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt12.priors',$
           prefix=path + 'fitresults/KELT-12.',$
           maxsteps=maxsteps,nthin=nthin,circular=[1],$
           debug=debug, verbose=verbose,nthread=nthread,$
           mistsedfile=path+'kelt12.sed',/fitslope,/ttv

end
