pro fitgj9827,debug=debug

maxsteps=10000
nthin=20
exofastv2, nplanets=3, tranpath='gj9827.Kepler.dat',fluxfile='gj9827.flux',$
           priorfile='gj9827.priors',debug=debug, prefix='gj9827.',maxsteps=maxsteps,$
           nthin=nthin,fittran=[1,1,1],fitrv=[0,0,0],circular=[1,1,1],/longcadence,/noyy,/noclaret,/earth

;stop

end
