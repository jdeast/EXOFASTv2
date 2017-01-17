;; returns the fitted parameters (do I care?)
function str2pars, str, scale=scale, name=name

  tofit = *(str.tofit)

  ;; slow way?
  pars = dblarr(n_elements(tofit[0,*]))
  for i=0, n_elements(pars)-1 do $
     pars[i] = str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value
     
  if arg_present(scale) then begin
     scale = dblarr(n_elements(tofit[0,*]))
     for i=0, n_elements(scale)-1 do $
        scale[i] = str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).scale
  endif

  if arg_present(name) then begin
     name = strarr(n_elements(tofit[0,*]))
     for i=0, n_elements(name)-1 do begin
        names = tag_names(str.(tofit[0,i])[tofit[1,i]])
        name[i] = names[tofit[2,i]] 
        if tofit[1,i] ne 0 then name[i] += '_' + strtrim(tofit[1,i],2)
     endfor
  endif

  return, pars

end
