;+
; NAME:
;   mkprior2
; PURPOSE:
;
;   Prints the contents of an EXOFASTv2 prior file so the next fit
;   starts at the most likely model from a previous EXOFASTv2
;   run. Often times such iteration can be helpful to make the fit
;   converge faster, since it is often difficult to find a good best
;   fit of complex fits with strong covariances and many parameters.
;
;   Note: If priors widths or bounds were supplied in the original
;   fit, these will be preserved and override the best-fit values.
;
; INPUTS:
;
;   FILENAME - The name of an IDL save file containing an MCMCSS
;              structure output by EXOFASTv2. Takes precendence over
;              MCMCSS.
;
;              Hint: a poorly behaving run can be terminated and still
;              generate this file by typing "ctrl+c", then
;              "!stopnow=1" and ".con" in the window.
;
;   MCMCSS   - The stellar stucture output by EXOFASTv2.
;
;   PRIORFILENAME - If specified, create a new priorfile with this
;                   name. Otherwise, output to the screen.
;
; REVISION HISTORY:
;   2018/03 - Public Release - Jason Eastman (CfA)
;
;-

function getpriorline2, parameter, ndx, num=num, ss=ss

label = parameter.label
if n_elements(num) ne 0 then begin
   if num ne 0 then label = label + '_' + strtrim(num,2)
endif

nlines = n_elements(*parameter.prior_new)
if nlines eq 0 then begin
   if parameter.fit then begin
      line = label + ' ' + string(parameter.best,format='(f0.30)')
      return, line
   endif else return, ''
endif

for i=0L, nlines-1 do begin

   prior = (*parameter.prior_new)[i]
 
   ;; if a map to a detrending variable
   if prior.value[4] ne -1 then begin
      priorval = (*(ss.(prior.value[0])))[prior.value[1]].(prior.value[2])[prior.value[3]].(prior.value[4])[prior.value[5]].label
   ;; if a map to another variable
   endif else if prior.value[3] ne -1 then begin
      priorval = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].label + '_' + strtrim(long(prior.value[1]),2)
   endif else if prior.value[2] ne -1 then begin
      priorval = ss.(prior.value[0])[prior.value[1]].label      
   endif else begin 
      priorval = strtrim(string(prior.value[0],format='(f0.30)'),2)
   endelse
   
   upperbound = prior.upperbound
   lowerbound = prior.lowerbound
   priorwidth = prior.gaussian_width
   bestval = strtrim(string(parameter.value[ndx],format='(f0.30)'),2)

   ;; if the parameter was fixed before, keep it fixed
   if priorwidth eq 0 then begin
      line = label + ' ' + priorval + ' 0'
   endif else if finite(priorwidth) and priorwidth gt 0d0 then begin
      ;; if a gaussian prior was supplied before, keep it 
      ;; (but still start at best fit)
      
      width = strtrim(priorwidth,2)
      if finite(upperbound) then upperbound = strtrim(upperbound,2) $
      else upperbound = 'Inf'
      if finite(lowerbound) then lowerbound = strtrim(lowerbound,2) $
      else lowerbound = '-Inf'
      line = label + ' ' + priorval + ' ' + width + ' ' + lowerbound + ' ' + upperbound + ' ' + bestval
      
   endif else if finite(upperbound) or finite(lowerbound) then begin
      ;; if just bounds were supplied, keep the bounds, but adjust the
      ;; starting value to the best fit value
      width = '-1'
      if finite(upperbound) then upperbound = strtrim(upperbound,2) $
      else upperbound = ''
      if finite(lowerbound) then lowerbound = strtrim(lowerbound,2) $
      else begin
         if upperbound eq '' then lowerbound = '' $
         else lowerbound = '-Inf'
      endelse
      line = label + ' ' + bestval + ' ' + width + ' ' + lowerbound + ' ' + upperbound
   endif else if parameter.fit then begin
      ;; otherwise, if it's a fitted parameter, start at the best value
      line = label + ' ' + bestval
   endif else line = '' ;; otherwise, don't bother

;print, label
;stop
endfor   

return, line

end

pro mkprior2, filename=filename, mcmcss=mcmcss, priorfilename=priorfilename

nargs = 0
if n_elements(filename) ne 0 then nargs++
if n_elements(mcmcss) ne 0 then nargs++
if n_elements(priorfilename) ne 0 then nargs++

;; for use without a license
if nargs lt 2 and (lmgr(/vm) or lmgr(/runtime)) then begin
   par = command_line_args(count=numargs)
   
   if numargs ne 2 then message, 'Must specify FILENAME and PRIORFILENAME'

   for i=0L, numargs-1 do begin
      if strpos(par[i],'=') ne -1 then begin
         entries = strsplit(par[i],'=',/extract)
         if strupcase(entries[0]) eq 'FILENAME' then filename = entries[1]
         if strupcase(entries[0]) eq 'PRIORFILENAME' then priorfilename = entries[1]
      endif
   endfor

endif


if n_elements(priorfilename) eq 1 then openw, lun, priorfilename, /get_lun $
else lun = -1   

if n_elements(filename) ne 0 then begin
   if ~file_test(filename) then message, filename + ' does not exist'
   restore, filename
endif else if n_elements(mcmcss) eq 0 then $
   message, 'Must specify either FILENAME or MCMCSS'

;; use the best model as the starting values
minchi2 = min(*mcmcss.chi2,ndx)

;; star
for i=0L, mcmcss.nstars-1 do begin
   printf, lun, '# star ' + mcmcss.star[i].label
   for j=0, n_tags(mcmcss.star[i])-1 do begin
      if (size(mcmcss.star[i].(j)))[2] eq 8 then begin
         line = getpriorline2(mcmcss.star[i].(j), ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; telescopes
for i=0L, mcmcss.ntel-1 do begin
   printf, lun, '# ' + mcmcss.telescope[i].label
   for j=0, n_tags(mcmcss.telescope[i])-1 do begin
      if (size(mcmcss.telescope[i].(j)))[2] eq 8 then begin        
         line = getpriorline2(mcmcss.telescope[i].(j), ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line 
      endif
   endfor

   if tag_exist((*(mcmcss.telescope[i].rvptrs)), 'NADD') then begin
      for j=0, (*(mcmcss.telescope[i].rvptrs)).nadd-1 do begin
         line = getpriorline2((*(mcmcss.telescope[i].rvptrs)).detrendaddpars[j],ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endfor
   endif

   if tag_exist((*(mcmcss.telescope[i].rvptrs)), 'NMULT') then begin
      for j=0, (*(mcmcss.telescope[i].rvptrs)).nmult-1 do begin
         line = getpriorline2((*(mcmcss.telescope[i].rvptrs)).detrendmultpars[j],ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endfor
   endif

endfor

;; astrometry
if tag_exist(mcmcss,'nastrom') then begin
   for i=0L, mcmcss.nastrom-1 do begin
      printf, lun, '# ' + mcmcss.astrom[i].label
      for j=0, n_tags(mcmcss.astrom[i])-1 do begin
         if (size(mcmcss.astrom[i].(j)))[2] eq 8 then begin        
            line = getpriorline2(mcmcss.astrom[i].(j), ndx, num=i, ss=mcmcss)
            if line ne '' then printf, lun, line 
         endif
      endfor
   endfor
endif

;; planets
for i=0L, mcmcss.nplanets-1 do begin
   printf, lun, '# ' + mcmcss.planet[i].label
   for j=0, n_tags(mcmcss.planet[i])-1 do begin
      if (size(mcmcss.planet[i].(j)))[2] eq 8 then begin         
         line = getpriorline2(mcmcss.planet[i].(j), ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; bands
for i=0L, mcmcss.nband-1 do begin
   printf, lun, '# ' + mcmcss.band[i].label
   for j=0, n_tags(mcmcss.band[i])-1 do begin
      if (size(mcmcss.band[i].(j)))[2] eq 8 then begin        
         line = getpriorline2(mcmcss.band[i].(j), ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; transits
for i=0L, mcmcss.ntran-1 do begin
   printf, lun, '# ' + mcmcss.transit[i].label
   for j=0, n_tags(mcmcss.transit[i])-1 do begin
      if (size(mcmcss.transit[i].(j)))[2] eq 8 then begin
         line = getpriorline2(mcmcss.transit[i].(j), ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endif
   endfor

   if tag_exist((*(mcmcss.transit[i].transitptrs)), 'NADD') then begin
      for j=0, (*(mcmcss.transit[i].transitptrs)).nadd-1 do begin
         line = getpriorline2((*(mcmcss.transit[i].transitptrs)).detrendaddpars[j],ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endfor
   endif

   if tag_exist((*(mcmcss.transit[i].transitptrs)), 'NMULT') then begin
      for j=0, (*(mcmcss.transit[i].transitptrs)).nmult-1 do begin
         line = getpriorline2((*(mcmcss.transit[i].transitptrs)).detrendmultpars[j],ndx, num=i, ss=mcmcss)
         if line ne '' then printf, lun, line
      endfor
   endif
   
endfor      

if lun ne -1 then free_lun, lun

end
