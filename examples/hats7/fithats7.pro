pro fithats7, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hats7'])

exofastv2, nplanets=1, rvpath=path+'HATS-7b.*.rv',tranpath=path+'n20*.dat',$
           priorfile=path+'hats7.priors',$
           prefix=path+'fitresults' + path_sep()+'HATS-7b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], exptime=[2,2],ninterp=[1,1],$
           mistsedfile='hats7.sed',$
           debug=debug,verbose=verbose, nthread=nthread,/nomist, /parsec

end
