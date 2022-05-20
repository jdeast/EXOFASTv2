pro fitkelt14,debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt14'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','kelt14','fitresults'])

exofastv2, nplanets=1, rvpath=path+'KELT-14b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt14.priors',$
           prefix=outpath+'KELT-14b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,nthread=nthread,$
           /nomist, /torres

end
