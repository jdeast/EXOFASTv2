;; This example showcases the BEER phase curve terms (BEaming,
;; Ellipsoidal, and Reflection) and uses a TESS LC
pro fittoi2165, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','toi2165'])

exofastv2, nplanets=1, tranpath=path + 'n20*.without_systematics.dat',$
           priorfile=path + 'toi2165.priors.3', $
           prefix=path + 'fitresults/TOI-2165b.', $
           maxsteps=maxsteps, $
           nthin=nthin, circular=[1], $
           debug=debug, verbose=verbose, $ 
           ntemps=8, mistsedfile = path + 'toi2165.sed',$
           fitreflect=['TESS'],fitthermal=['TESS'],$
           fitellip=['TESS'], fitbeam=[1]
           exptime=[30],ninterp=[10]

end
