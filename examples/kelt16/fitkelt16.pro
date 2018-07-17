pro fitkelt16, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps=5000
if n_elements(nthin) eq 0 then nthin=10

exofastv2, nplanets=1, rvpath='KELT-16b.*.rv',tranpath='n20*.dat',$
           priorfile='kelt16.priors',$
           prefix='FitResults/KELT-16b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres;,fluxfile='kelt16.sed'

end
