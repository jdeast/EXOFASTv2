pro fithat3, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
             nthread=nthread, outpath=outpath

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

if n_elements(outpath) eq 0 then $
   outpath = 'modeling' + path_sep() + 'hat3' + path_sep() + 'fitresults' + path_sep()

;; Fit using the Torres relation
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', $
           rvpath=path+'HAT-3b.*.rv',priorfile=path+'hat3.torres.priors', $
           prefix=outpath +'HAT-3b.Torres.', maxsteps=maxsteps, $
           nthin=nthin, circular=[1], fitrv=[1],fittran=[1], $
           debug=debug, verbose=verbose, /nomist, /torres, nthread=nthread


end
