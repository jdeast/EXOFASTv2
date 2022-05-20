pro fitgj9827, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','gj9827'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','gj9827','fitresults'])

exofastv2, nplanets=3, tranpath=path+'n20090101.Kepler.K2.dat',$
           mistsedfile=path+'gj9827.flux',$
           priorfile=path+'gj9827.priors',debug=debug, verbose=verbose, $
           prefix=outpath+'GJ9827.',maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1,1],fitrv=[0,0,0],circular=[1,1,1],$
           /longcadence,/noclaret,/earth,/nomist,nthread=nthread

end
