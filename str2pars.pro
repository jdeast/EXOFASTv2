;; returns the fitted parameters (do I care?)
function str2pars, str, scale=scale, name=name, angular=angular,ndx=ndx

  if n_elements(ndx) eq 0 then ndx=0L
  angular = [-1]
  tofit = *(str.tofit)

  ;; slow way?
  pars = dblarr(n_elements(tofit[0,*]))
  for i=0, n_elements(pars)-1 do begin
     if tofit[3,i] eq -1 then begin
        pars[i] = str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value[ndx]
        if str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).unit eq 'Radians' then angular = [angular,i]
     endif else begin
        pars[i] = (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value[ndx]
        if (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].unit eq 'Radians' then angular = [angular,i]
     endelse
  endfor

  if arg_present(scale) then begin
     scale = dblarr(n_elements(tofit[0,*]))
     for i=0, n_elements(scale)-1 do begin
        if tofit[3,i] eq -1 then begin
           scale[i] = str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).scale
        endif else begin
           scale[i] = (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].scale
        endelse
     endfor
  endif

  if arg_present(name) then begin
     name = strarr(n_elements(tofit[0,*]))
     for i=0, n_elements(name)-1 do begin
        if tofit[3,i] eq -1 then begin
           names = tag_names(str.(tofit[0,i])[tofit[1,i]])
           name[i] = names[tofit[2,i]] 
        endif else begin
           name[i] = (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].label
        endelse

        if tofit[1,i] ne 0 then name[i] += '_' + strtrim(tofit[1,i],2)
     endfor
  endif

  if n_elements(angular) gt 1 then angular = angular[1:*]
  return, pars

end
