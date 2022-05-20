pro fitkelt6, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt6'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','kelt6','fitresults'])

exofastv2, nplanets=2, rvpath=path+'KELT-6b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'kelt6.priors',$
           prefix=outpath+'KELT-6b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1,1],fittran=[1,0],circular=[0,0], $
           debug=debug,verbose=verbose,$
           mistsedfile=path+'kelt6.sed', nthread=nthread

end
