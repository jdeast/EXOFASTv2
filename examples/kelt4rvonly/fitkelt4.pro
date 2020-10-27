pro fitkelt4, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','kelt4rvonly'])

;; NOTE: kelt4.priors.2 was created using mkprior on an initial run
;; using kelt4.priors.

exofastv2, nplanets=1, rvpath=path+'KELT-4b.*.rv',$
           priorfile=path+'kelt4.priors.2',$
           prefix=path+'fitresults' + path_sep()+'KELT-4Ab.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],circular=[0], $
           debug=debug,verbose=verbose, nthread=nthread

end
