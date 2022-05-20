;; This example compares several different methods of determining the
;; stellar properties

pro fithat3_comparestar, debug=debug, verbose=verbose, dontstop=dontstop, maxsteps=maxsteps, nthin=nthin, outpath=outpath

if n_elements(nthin) eq 0 then nthin=1
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

if n_elements(outpath) eq 0 then $
   outpath = filepath('',root_dir=getenv('HOME'),subdir=['modeling','hat3','fitresults'])

;; Fit using SED+MIST with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.SED+MIST+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, mistsedfile=path+'hat3.sed', dontstop=dontstop

;; Fit using the Torres relation with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.Torres+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, /nomist, /torres, dontstop=dontstop

;; Fit using YY isochrones with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.YY+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, dontstop=dontstop, /yy, /nomist

;; Fit using MIST isochrones with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.MIST+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, dontstop=dontstop, mistsedfile=path+'hat3.sed'

;; Fit using SED with the transit
exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', rvpath=path+'HAT-3b.HIRES.rv',priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.SED+Transit.', maxsteps=maxsteps, nthin=nthin, circular=[1], fitrv=[1],fittran=[1], debug=debug, verbose=verbose, mistsedfile=path+'hat3.sed', /nomist, dontstop=dontstop

;; Fit only the star with MIST and the SED (no transit or RV)
exofastv2, nplanets=0, priorfile=path+'hat3.priors', prefix=outpath + 'HAT-3b.SED+MIST.', maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose, mistsedfile=path+'hat3.sed', dontstop=dontstop

comparestars, outpath + 'HAT-3b.*.mcmc.idl', psname=path+'comparestars.ps'

end
