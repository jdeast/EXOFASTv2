function mcmcss2ss, mcmcss, filename=filename

if n_elements(filename) eq 1 then restore, filename

;; backward compatibility                              
if tag_exist(mcmcss,'requiresecondary') then begin
   requiresecondary = mcmcss.requiresecondary
endif else begin
   requiresecondary = 0
endelse

if tag_exist(mcmcss,'noprimary') then begin
   noprimary = mcmcss.noprimary
endif else begin
   noprimary = 0
endelse

if tag_exist(mcmcss.transit,'claret') then begin
   claret = mcmcss.transit.claret
endif else if tag_exist(mcmcss.band,'claret') then begin
   claret = mcmcss.band.claret
endif else if tag_exist(mcmcss,'claret') then begin
   claret = mcmcss.claret
endif else begin
   claret = 0
endelse

if tag_exist(mcmcss,'seddeblend') then begin
   seddeblend = mcmcss.seddeblend
endif else if tag_exist(mcmcss,'diluted') then begin
   seddeblend = mcmcss.diluted
endif else begin
   seddeblend = 0
endelse


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
          noprimary=noprimary, requiresecondary=requiresecondary,$
          ninterp=mcmcss.ninterp, exptime=mcmcss.exptime, earth=mcmcss.earth, i180=mcmcss.i180,$
          seddeblend=seddeblend, fitdilute=fitdilute, $
          fitthermal=mcmcss.fitthermal, fitreflect=mcmcss.fitreflect, fitphase=mcmcss.fitphase, $
          fitellip=mcmcss.fitellip, fitbeam=mcmcss.fitbeam, derivebeam=mcmcss.derivebeam, $
          chen=mcmcss.chen, yy=mcmcss.yy, torres=mcmcss.torres, nomist=~mcmcss.mist, parsec=mcmcss.parsec, $
          noclaret=~claret, alloworbitcrossing=alloworbitcrossing, logname=mcmcss.logname,$
          fitspline=mcmcss.fitspline, splinespace=mcmcss.splinespace, fitwavelet=mcmcss.fitwavelet, $
          novcve=mcmcss.novcve, nochord=mcmcss.nochord, fitsign=mcmcss.fitsign, $
          chi2func=mcmcss.chi2func, fittt=mcmcss.fittt, rvepoch=mcmcss.rvepoch,delay=mcmcss.delay,prefix=mcmcss.prefix, $
          sedrange=mcmcss.sedrange, transitrange=mcmcss.transitrange, rvrange=mcmcss.rvrange, emrange=mcmcss.emrange)

return, ss

end
