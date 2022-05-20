pro fithat3notran, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, outpath=outpath

if n_elements(nthin) eq 0 then nthin=1
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','hat3','fitresults'])


;; Fit using the Torres relation
exofastv2, nplanets=1, $;tranpath=path+'n20070428.Sloani.KepCam.dat', $
           rvpath=path+'HAT-3b.*.rv',priorfile=path+'hat3.torres.priors', $
           prefix=outpath+'HAT-3b.notran.', maxsteps=maxsteps, $
           nthin=nthin, circular=[1], fitrv=[1],$;fittran=[1], $
           debug=debug, verbose=verbose, /nomist, /torres

end
