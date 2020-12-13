pro printdeltavdiff, idlfile

restore, idlfile
chi2 = reform((*mcmcss.chi2),mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains)

for i=0L, mcmcss.ntran-1L do begin
   for j=i+1L, mcmcss.ntran-1L do begin
      diff = (reform(mcmcss.transit[i].tdeltav.value-mcmcss.transit[j].tdeltav.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]
      summarizepar, diff, value=value, errlo=errlo, errhi=errhi
      if value gt 0d0 then sigdiff = double(value)/double(errlo) $
      else sigdiff = double(value)/double(errhi)
      print, mcmcss.transit[i].label + '-' + mcmcss.transit[j].label + ': ' + strtrim(value,2) + ' +' + strtrim(errhi,2) + ' /-' + strtrim(errlo,2) + ' (' + strtrim(sigdiff,2) + ' sigma significant)'
   endfor
endfor

end
