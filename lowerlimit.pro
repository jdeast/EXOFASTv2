pro lowerlimit, filename, names=names, sigma=sigma

if n_elements(sigma) eq 0 then sigma = 3d0
if n_elements(names) eq 0 then $
   names = ['R_P','R_P/R_*','\rho_P','logg_P','i','b']

sigmapercent = erf(sigma/sqrt(2d0))
restore, filename

chi2 = reform((*mcmcss.chi2),mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains)

for n=0, n_elements(names)-1 do begin
   found = 0

   for i=0, n_tags(mcmcss)-1 do begin
      for j=0, n_elements(mcmcss.(i))-1 do begin
         for k=0, n_tags(mcmcss.(i)[j])-1 do begin
            
            ;; this captures the detrending variables
            if (size(mcmcss.(i)[j].(k)))[1] eq 10 then begin
               if ptr_valid(mcmcss.(i)[j].(k)) then begin
                  for l=0L, n_tags(*(mcmcss.(i)[j].(k)))-1 do begin
                     if (size((*(mcmcss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                        for m=0L, n_elements((*(mcmcss.(i)[j].(k))).(l))-1 do begin
                           if tag_exist((*(mcmcss.(i)[j].(k))).(l)[m],'latex') then begin
                              if (*(mcmcss.(i)[j].(k))).(l)[m].latex eq names[n] then begin
                                 found = 1
                                 units = (*(mcmcss.(i)[j].(k))).(l)[m].unit
                                 par = (reform((*(mcmcss.(i)[j].(k))).(l)[m].value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]

                                 nan = where(finite(par,/nan))
                                 if nan[0] ne -1 then par[nan] = !values.d_infinity

                                 nsteps = n_elements(par)
                                 lowerlimit = (par[reverse(sort(par))])[round(sigmapercent*nsteps)]
                                 print, strtrim(sigma,2) + ' sigma lower limit on ' + names[n] + ',' + strtrim(j,2) + ' (' + units + ') ' + strtrim(lowerlimit,2)
                                 found = 1
                              endif
                           endif
                        endfor
                     endif
                  endfor
               endif
            ;; this captures everything else
            endif else if n_tags(mcmcss.(i)[j].(k)) ne 0 then begin
               if n_tags(mcmcss.(i)[j].(k)) ne 0 then begin
                  if tag_exist(mcmcss.(i)[j].(k),'latex') then begin
                     if mcmcss.(i)[j].(k).latex eq names[n] then begin
                        found = 1
                        units = mcmcss.(i)[j].(k).unit
                        if units eq '' then unitstr = ' ' $
                        else unitstr = ' (' + units + ') '

                        par = (reform(mcmcss.(i)[j].(k).value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains]

                        nan = where(finite(par,/nan))
                        if nan[0] ne -1 then par[nan] = !values.d_infinity

                        nsteps = n_elements(par)
                        lowerlimit = (par[reverse(sort(par))])[round(sigmapercent*nsteps)]
                        print, strtrim(sigma,2) + ' sigma lower limit on ' + names[n] + ',' + strtrim(j,2) + unitstr + strtrim(lowerlimit,2)      
                        found = 1
                     endif
                  endif
               endif
            endif      
         endfor
      endfor
   endfor
   
   if ~found then print, 'no match found for ' + names[n]
  
endfor

end
