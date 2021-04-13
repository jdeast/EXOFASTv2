pro fit8361, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','ep247098361'])

if n_elements(outpath) eq 0 then $
   outpath = 'modeling' + path_sep() + 'ep8361' + path_sep() + 'fitresults' + path_sep()

exofastv2, nplanets=1, tranpath=path+'ep247098361.Kepler.dat',rvpath=path+'ep247098361.APF.rv',$
           priorfile=path+'ep247098361.priors',debug=debug, verbose=verbose, $
           prefix=outpath + 'EPIC247098361b.',$
           maxsteps=maxsteps,nthin=nthin,$
           fittran=[1],fitrv=[1],circular=[0],/longcadence, nthread=nthread

end
