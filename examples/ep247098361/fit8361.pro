pro fit8361, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps=1000
if n_elements(nthin) eq 0 then nthin=1

exofastv2, nplanets=1, tranpath='ep247098361.Kepler.dat',rvpath='ep247098361.APF.rv',$
           priorfile='ep247098361_priors.dat',debug=debug, verbose=verbose, $
           prefix='FitResults/ep247098361.MIST.',$
           maxsteps=maxsteps,nthin=nthin,$
           fittran=[1],fitrv=[1],circular=[0],/longcadence

end
