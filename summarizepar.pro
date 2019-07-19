pro summarizepar, pars, label=label, unit=unit, latex=latex, best=best, logname=logname, value=value, errlo=errlo, errhi=errhi, medianpars=medianpars, charsize=charsize, mode=mode, noplot=noplot

if n_elements(unit) eq 0 then unit = ''
if n_elements(label) eq 0 then label = ''
if n_elements(latex) eq 0 then latex = ''
if n_elements(best) eq 0 then best = !values.d_infinity

sz = size(pars)
nchains = sz[2]

;; check for bad values
;bad = where(~finite(pars),complement=good)
good = where(finite(pars),complement=bad,nsteps)
if bad[0] ne -1 then printandlog, $
   "Warning: NaNs in " + label + " distribution; NaNs will be ignored in all summaries",logname

if nsteps eq 0 then begin
   printandlog,"ERROR: all NaNs in " + label + " distribution; skipping parameter",logname
   return
endif

;; if angular, center distribution about the mode first
if strupcase(unit) eq 'DEGREES' then halfrange = 180d0 $
else if strupcase(unit) eq 'RADIANS' then halfrange = !dpi
if strupcase(unit) eq 'DEGREES' or strupcase(unit) eq 'RADIANS' then begin
   ;; find the mode
   hist = histogram(pars[good],nbins=100,locations=x,/nan)
   max = max(hist,modendx)
   modepar = x[modendx]
   
   toohigh = where(pars[good] gt (modepar + halfrange))
   if toohigh[0] ne -1 then pars[good[toohigh]] -= 2.d0*halfrange
   
   toolow = where(pars[good] lt (modepar - halfrange))
   if toolow[0] ne -1 then pars[good[toolow]] += 2.d0*halfrange
endif                     

sorted = sort(pars[good])
halfsigma = erf(1d/sqrt(2d0))/2d

if keyword_set(mode) then begin
;   mode = kde(
;   sigma = 2d0*halfsigma
;
;   repeat begin
;      maxprob =    
;   endrep until 1
endif else begin
   ;; median and 68% confidence interval
   medndx = nsteps/2d0
   lowsigndx = round(nsteps/2.d0 - nsteps*halfsigma)
   hisigndx = round(nsteps/2.d0 + nsteps*halfsigma)
   

   medvalue = pars[good[sorted[medndx]]]
   upper = pars[good[sorted[hisigndx]]] - medvalue
   lower = medvalue - pars[good[sorted[lowsigndx]]]
endelse

if n_elements(medianpars) eq 1 then medianpars = [medvalue,upper,lower] $
else medianpars = [[medianpars],[[medvalue,upper,lower]]]

xmax = (medvalue + 4*upper) < max(pars[good])
xmin = (medvalue - 4*lower) > min(pars[good])

if xmin eq xmax then begin
   printandlog, 'WARNING: ' + label + ' is singularly valued.',logname
endif else if ~keyword_set(noplot) then begin
               
   ;; plot labels
   xtitle='!3' + exofast_textoidl(latex)
   ytitle='!3Probability'

   hist = histogram(pars[good],nbins=100,locations=x,min=xmin,max=xmax,/nan)
   plot, x, hist/double(total(hist)), psym=10, xtitle=xtitle, ytitle=ytitle,$
         charsize=charsize,xstyle=1, xrange=[xmin,xmax], font=1
   
   for l=0L, nchains-1L do  begin
      hist = histogram(pars[*,l],nbins=100,locations=x,min=xmin,max=xmax,/nan)
      if total(hist) gt 0 then $
         oplot, x, hist/double(total(hist)), psym=10,color=l*255d0/nchains
   endfor

   hist = histogram(pars[good],nbins=100,locations=x,min=xmin,max=xmax,/nan)
   oplot, x, hist/double(total(hist)), psym=10, thick=3
   
   ;; if the best parameters are given, overplot them on the PDFs
   if finite(best) then $
      oplot, [best,best],[-9d9,9d9]

endif

;; format values for table (rounded appropriately)
;; round the high error to 2 sig figs
exphi=fix(alog10(upper))
if (upper lt 1d0) then exphi=exphi-1
roundhi=round(upper/10.d0^(exphi-1d0),/L64)*10.d0^(exphi-1d0)
if (roundhi gt 10) then errhi = strtrim(round(roundhi,/L64),2) $
else errhi = string(roundhi,format='(f255.'+strtrim(1-exphi,2)+')')

;; round the low error to 2 sig figs
explo=fix(alog10(lower))
if (lower lt 1d0) then explo=explo-1
roundlo=round(lower/10.d0^(explo-1d0),/L64)*10.d0^(explo-1d0)
if (roundlo gt 10) then errlo = strtrim(round(roundlo,/L64),2) $
else errlo = string(roundlo,format='(f255.'+strtrim(1-explo,2)+')')

;; round the value to the greater number of sig figs
ndec = long(1 - (exphi < explo))
if ndec eq 0 then value = string(medvalue,format='(i255)') $
else if ndec lt 0 then $
   value=round(round(medvalue/10.d0^(-ndec),/L64)*10.d0^(-ndec),/L64) $
else value = string(medvalue,format='(f255.'+strtrim(ndec,2)+')')

end
