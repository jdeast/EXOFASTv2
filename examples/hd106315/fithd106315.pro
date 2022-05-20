pro fithd106315, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hd106315'])

if n_elements(outpath) eq 0 then $
      outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','hd106315','fitresults'])

exofastv2, nplanets=2, tranpath=path + 'n20090101.Kepler.K2.dat',mistsedfile=path+'hd106315.flux',$
           priorfile=path + 'hd106315.priors',debug=debug, verbose=verbose, $
           prefix=outpath + 'hd106315.',maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1],fitrv=[0,0],circular=[0,0],/longcadence, nthread=nthread

;stop

end
