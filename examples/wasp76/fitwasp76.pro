pro fitwasp76, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps=5000
if n_elements(nthin) eq 0 then nthin=10

exofastv2, nplanets=1, rvpath='WASP-76b.*.rv',tranpath='n20*.dat',$
           priorfile='wasp76.priors',$
           prefix='FitResults/WASP-76b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres;,fluxfile='wasp76.sed'

end
