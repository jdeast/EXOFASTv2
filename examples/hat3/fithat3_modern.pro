pro fithat3_modern, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
             nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','hat3','fitresults'])

exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', $
           rvpath=path+'HAT-3b.*.rv',$
           priorfile=path+'HAT-3.priors', $
           mistsedfile=path+'HAT-3.sed',$
           prefix=outpath +'HAT-3b.', $
           maxsteps=maxsteps, nthin=nthin, $
           debug=debug, verbose=verbose, $
           nthread=nthread, $
           fehsedfloor=0.08d0

end
