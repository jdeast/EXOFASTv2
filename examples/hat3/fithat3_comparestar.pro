;; This example compares several different methods of determining the
;; stellar properties

pro fithat3_comparestar, debug=debug, verbose=verbose, dontstop=dontstop, maxsteps=maxsteps, nthin=nthin

;; default to a very short run (not well-mixed or publication quality)
if n_elements(maxsteps) eq 0 then maxsteps=100 ;; 50000
if n_elements(nthin) eq 0 then nthin=1 ;; 50
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

;; Fit using the Torres relation with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=path+'HAT-3b.Torres+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, /noyy, /torres, dontstop=dontstop

;; Fit using YY isochrones with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=path+'HAT-3b.YY+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, dontstop=dontstop

;; Fit using MIST isochrones with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=path+'HAT-3b.MIST+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, /noyy, /mist, dontstop=dontstop, fluxfile=path+'hat3.flux'

;; Fit using SED (current precision) with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=path+'HAT-3b.SED+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, dontstop=dontstop

;; Fit using SED+MIST (current precision) with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=path+'HAT-3b.SED+MIST+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, /mist, dontstop=dontstop

;; Fit using SED with 5 uas precision with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.gaia.priors', prefix=path+'HAT-3b.SED+5uas+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, dontstop=dontstop

;; Fit using SED with 5 uas precision + MIST with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.gaia.priors', prefix=path+'HAT-3b.SED+5uas+MIST+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, /mist, dontstop=dontstop

;; Fit only the star with MIST and the SED (no transit or RV)
exofastv2, nplanets=0, priorfile=path+'hat3.priors', prefix=path+'HAT-3b.SED+MIST.', maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, /mist, dontstop=dontstop

;; Fit only the star with MIST and the SED with 5 uas precision (no transit or RV)
exofastv2, nplanets=0, priorfile=path+'hat3.gaia.priors', prefix=path+'HAT-3b.SED+5uas+MIST.', maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose, fluxfile=path+'hat3.flux', /noyy, /mist, dontstop=dontstop

comparestars, path+'HAT-3b.*.mcmc.idl', psname=path+'comparestars.ps'

end
