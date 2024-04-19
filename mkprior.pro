;+
; NAME:
;   mkprior
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

function getpriorline, parameter, ndx, num=num

label = parameter.label
if n_elements(num) ne 0 then begin
   if num ne 0 then label = label + '_' + strtrim(num,2)
endif

bestval = strtrim(string(parameter.value[ndx],format='(f0.30)'),2)

;; if the parameter was fixed before, keep it fixed
if parameter.priorwidth eq 0 then begin
   priorval = strtrim(string(parameter.prior,format='(f0.30)'),2)
   line = label + ' ' + priorval + ' 0'
endif else if finite(parameter.priorwidth) then begin
   ;; if a gaussian prior was supplied before, keep it 
   ;; (but still start at best fit)
   priorval = strtrim(string(parameter.prior,format='(f0.30)'),2)

   width = strtrim(parameter.priorwidth,2)
   if finite(parameter.upperbound) then upperbound = strtrim(parameter.upperbound,2) $
   else upperbound = 'Inf'
   if finite(parameter.lowerbound) then lowerbound = strtrim(parameter.lowerbound,2) $
   else lowerbound = '-Inf'

   if double(width) lt 0d0 then width = '-1'

   if width eq '-1' and upperbound eq 'Inf' and lowerbound eq '-Inf' then begin
      line = label + ' ' + bestval
   endif else if width eq '-1' then begin
      line = label + ' ' + bestval + ' ' + width + ' ' + lowerbound + ' ' + upperbound
   endif else begin
      line = label + ' ' + priorval + ' ' + width + ' ' + lowerbound + ' ' + upperbound + ' ' + bestval
   endelse
   
endif else if finite(parameter.upperbound) or finite(parameter.lowerbound) then begin
   ;; if just bounds were supplied, keep the bounds, but adjust the
   ;; starting value to the best fit value
   width = '-1'
   if finite(parameter.upperbound) then upperbound = strtrim(parameter.upperbound,2) $
   else upperbound = ''
   if finite(parameter.lowerbound) then lowerbound = strtrim(parameter.lowerbound,2) $
   else begin
      if upperbound eq '' then lowerbound = '' $
      else lowerbound = '-Inf'
   endelse
   line = label + ' ' + bestval + ' ' + width + ' ' + lowerbound + ' ' + upperbound
endif else if parameter.fit then begin
   ;; otherwise, if it's a fitted parameter, start at the best value
   line = label + ' ' + bestval
endif else line = ''

return, line

end

pro mkprior, filename=filename, mcmcss=mcmcss, priorfilename=priorfilename

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
   for j=0, n_tags(mcmcss.star[i])-1 do begin
      if (size(mcmcss.star[i].(j)))[2] eq 8 then begin
         line = getpriorline(mcmcss.star[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; telescopes
for i=0L, mcmcss.ntel-1 do begin
   printf, lun, '# ' + mcmcss.telescope[i].label
   for j=0, n_tags(mcmcss.telescope[i])-1 do begin
      if (size(mcmcss.telescope[i].(j)))[2] eq 8 then begin        
         line = getpriorline(mcmcss.telescope[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line 
      endif
   endfor
endfor

;; astrometry
if tag_exist(mcmcss,'nastrom') then begin
   for i=0L, mcmcss.nastrom-1 do begin
      printf, lun, '# ' + mcmcss.astrom[i].label
      for j=0, n_tags(mcmcss.astrom[i])-1 do begin
         if (size(mcmcss.astrom[i].(j)))[2] eq 8 then begin        
            line = getpriorline(mcmcss.astrom[i].(j), ndx, num=i)
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
         line = getpriorline(mcmcss.planet[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; bands
for i=0L, mcmcss.nband-1 do begin
   printf, lun, '# ' + mcmcss.band[i].label
   for j=0, n_tags(mcmcss.band[i])-1 do begin
      if (size(mcmcss.band[i].(j)))[2] eq 8 then begin        
         line = getpriorline(mcmcss.band[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; spectrophotometry files
for i=0L, mcmcss.nspecfiles-1 do begin
   printf, lun, '# ' + mcmcss.specphot[i].label
   for j=0, n_tags(mcmcss.specphot[i])-1 do begin
      if (size(mcmcss.specphot[i].(j)))[2] eq 8 then begin        
         line = getpriorline(mcmcss.specphot[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line
      endif
   endfor
endfor

;; transits
for i=0L, mcmcss.ntran-1 do begin
   printf, lun, '# ' + mcmcss.transit[i].label
   for j=0, n_tags(mcmcss.transit[i])-1 do begin
      if (size(mcmcss.transit[i].(j)))[2] eq 8 then begin
         line = getpriorline(mcmcss.transit[i].(j), ndx, num=i)
         if line ne '' then printf, lun, line
      endif
   endfor

   if tag_exist((*(mcmcss.transit[i].transitptrs)), 'NADD') then begin
      for j=0, (*(mcmcss.transit[i].transitptrs)).nadd-1 do begin
         line = getpriorline((*(mcmcss.transit[i].transitptrs)).detrendaddpars[j],ndx, num=i)
         if line ne '' then printf, lun, line
      endfor
   endif

   if tag_exist((*(mcmcss.transit[i].transitptrs)), 'NMULT') then begin
      for j=0, (*(mcmcss.transit[i].transitptrs)).nmult-1 do begin
         line = getpriorline((*(mcmcss.transit[i].transitptrs)).detrendmultpars[j],ndx, num=i)
         if line ne '' then printf, lun, line
      endfor
   endif
   
endfor      

if lun ne -1 then free_lun, lun

end
