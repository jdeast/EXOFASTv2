function mcmcss2ss, mcmcss, filename=filename

if n_elements(filename) eq 1 then restore, filename

ss = mkss(rvpath=mcmcss.rvpath, tranpath=mcmcss.tranpath, astrompath=mcmcss.astrompath, dtpath=mcmcss.dtpath, $
          fluxfile=mcmcss.fluxfile, mistsedfile=mcmcss.mistsedfile,$
          fbolsedfloor=mcmcss.fbolsedfloor,teffsedfloor=mcmcss.teffsedfloor, fehsedfloor=mcmcss.fehsedfloor, $
          teffemfloor=mcmcss.teffemfloor, fehemfloor=mcmcss.fehemfloor, rstaremfloor=mcmcss.rstaremfloor,ageemfloor=mcmcss.ageemfloor,$
          oned=mcmcss.oned,nplanets=mcmcss.nplanets, nstars=mcmcss.nstars, $
          starndx=mcmcss.planet.starndx,$
          debug=mcmcss.debug, verbose=mcmcss.verbose, priorfile=mcmcss.priorfile, $
          fitrv=mcmcss.fitrv, fittran=mcmcss.fittran, fitdt=mcmcss.fitdt,rossiter=mcmcss.rossiter,$
          fitlogmp=mcmcss.fitlogmp,$
          circular=mcmcss.circular,fitslope=mcmcss.fitslope, fitquad=mcmcss.fitquad,tides=mcmcss.tides,$
          ttvs=mcmcss.ttvs, tivs=mcmcss.tivs, tdeltavs=mcmcss.tdeltavs,$
          longcadence=mcmcss.longcadence,rejectflatmodel=mcmcss.rejectflatmodel,$
          ninterp=mcmcss.ninterp, exptime=mcmcss.exptime, earth=mcmcss.earth, i180=mcmcss.i180,$
          diluted=mcmcss.diluted, fitdilute=mcmcss.fitdilute, $
          fitthermal=mcmcss.fitthermal, fitreflect=mcmcss.fitreflect, fitphase=mcmcss.fitphase, $
          fitellip=mcmcss.fitellip, fitbeam=mcmcss.fitbeam, derivebeam=mcmcss.derivebeam, $
          chen=mcmcss.chen, yy=mcmcss.yy, torres=mcmcss.torres, nomist=~mcmcss.mist, parsec=mcmcss.parsec, $
          noclaret=~mcmcss.claret, alloworbitcrossing=mcmcss.alloworbitcrossing, logname=mcmcss.logname,$
          fitspline=mcmcss.fitspline, splinespace=mcmcss.splinespace, fitwavelet=mcmcss.fitwavelet, $
          novcve=mcmcss.novcve, nochord=mcmcss.nochord, fitsign=mcmcss.fitsign, $
          chi2func=mcmcss.chi2func, fittt=mcmcss.fittt, rvepoch=mcmcss.rvepoch,delay=mcmcss.delay,prefix=mcmcss.prefix, $
          sedrange=mcmcss.sedrange, transitrange=mcmcss.transitrange, rvrange=mcmcss.rvrange, emrange=mcmcss.emrange)

return, ss

end
