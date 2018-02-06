pro fitkelt4, debug=debug

maxsteps=100
nthin=1

exofastv2, nplanets=1, rvpath='KELT-4b.*.rv',priorfile='kelt4.priors',prefix='KELT-4.',maxsteps=maxsteps,nthin=nthin,fitrv=[1],circular=[0], debug=debug

end
