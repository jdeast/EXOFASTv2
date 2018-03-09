;; This example compares several different methods of determining the
;; stellar properties

pro fithat3_comparestar, debug=debug, verbose=verbose, dontstop=dontstop

maxsteps=50000
nthin=5
dontstop=1 ;; don't stop when well mixed... makes them pretty, but takes a long time
;maxsteps=100
;nthin=1

;; Fit using the Torres relation with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-3b.Torres+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, /noyy, /torres, dontstop=dontstop

;; Fit using YY isochrones with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-3b.YY+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, dontstop=dontstop

;; Fit using MIST isochrones with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-3b.MIST+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, /noyy, /mist, dontstop=dontstop, fluxfile='hat3.flux'

;; Fit using SED (current precision) with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-3b.SED+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile='hat3.flux', /noyy, dontstop=dontstop

;; Fit using SED+MIST (current precision) with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.priors', prefix='HAT-3b.SED+MIST+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile='hat3.flux', /noyy, /mist, dontstop=dontstop

;; Fit using SED with 5 uas precision with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.gaia.priors', prefix='HAT-3b.SED+5uas+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile='hat3.flux', /noyy, dontstop=dontstop

;; Fit using SED with 5 uas precision + MIST with the transit
exofastv2, nplanets=1, tranpath='n20070428.Sloani.KepCam.dat', rvpath='hat3.rv',priorfile='hat3.gaia.priors', prefix='HAT-3b.SED+5uas+MIST+Transit.long.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile='hat3.flux', /noyy, /mist, dontstop=dontstop

;; Fit only the star with MIST and the SED (no transit or RV)
exofastv2, nplanets=0, priorfile='hat3.priors', prefix='HAT-3b.SED+MIST.long.', maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose, fluxfile='hat3.flux', /noyy, /mist, dontstop=dontstop

comparestars, 'HAT-3b.*.long.mcmc.idl', psname='comparestars.ps'

end
