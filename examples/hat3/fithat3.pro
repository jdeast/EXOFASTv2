pro fithat3, debug=debug

maxsteps=50000
nthin=2
;maxsteps=100
;nthin=1

exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-P-3b.torres4.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, /noyy, /torres

end
