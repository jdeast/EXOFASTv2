pro fithd106315,debug=debug

maxsteps=1d4
nthin=1d2
;maxsteps=1d5/3
;nthin=10
maxsteps = 100
nthin = 1

exofastv2, nplanets=2, tranpath='hd106315.Kepler.dat',fluxfile='hd106315.flux',$
           priorfile='hd106315.priors',debug=debug, prefix='hd106315.',$
           maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1],fitrv=[0,0],circular=[0,0],/longcadence

;stop

end
