pro fitwasp18, debug=debug, verbose=verbose,  maxsteps=maxsteps, nthin=nthin, nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','wasp18'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','wasp18','fitresults'])

exofastv2, nplanets=1, tranpath=path + 'n*.dat', rvpath=path+'*.rv',$ 
           priorfile=path+'wasp18.priors.2', mistsedfile=path+'wasp18.sed',$
           prefix=outpath + 'WASP-18b.', maxsteps=maxsteps, nthin=nthin, $
           debug=debug, verbose=verbose, ntemp=8,nthread=nthread,$
           fitthermal=['TESS'], fitreflect=['TESS'], fitellip=['TESS'],fitbeam=[1], fitphase=['TESS']
end
