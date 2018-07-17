pro fitkelt14,debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps=50000
if n_elements(nthin) eq 0 then nthin=50

exofastv2, nplanets=1, rvpath='KELT-14b.*.rv',tranpath='n20*.dat',$
           priorfile='kelt14.priors',$
           prefix='FitResults/KELT-14b.',maxsteps=maxsteps,$
           nthin=nthin,fitrv=[1],fittran=[1],circular=[0], $
           debug=debug,verbose=verbose,$
           /nomist, /torres;,fluxfile='kelt16.sed'

end
