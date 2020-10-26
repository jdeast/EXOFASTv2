pro fitkelt1, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt1'])

exofastv2, nplanets=1, rvpath=path+'KELT-1b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt1.priors',$
           prefix=path+'fitresults' + path_sep()+'KELT-1b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres;,mistsedfile=path+'kelt1.sed'

end
