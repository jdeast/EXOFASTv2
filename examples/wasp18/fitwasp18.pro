pro fitwasp18, debug=debug, verbose=verbose,  maxsteps=maxsteps, nthin=nthin, nthread=nthread

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','wasp18'])

exofastv2, nplanets=1, tranpath=path + 'n*.dat', rvpath=path+'*.rv',$ 
           priorfile=path+'wasp18.priors.2', mistsedfile=path+'wasp18.sed',$
           prefix=path+'fitresults/WASP-18b.', maxsteps=maxsteps, nthin=nthin, $
           debug=debug, verbose=verbose, ntemp=8,nthread=nthread,$
           fitthermal=['TESS'], fitreflect=['TESS'], fitbeam=[1]
end
