;; returns the fitted parameters (do I care?)
function str2parsarr, str, scale=scale, name=name, angular=angular

  angular = []
  tofit = *(str.tofit)

  ;; slow way?
  pars = dblarr(n_elements(tofit[0,*]),str.nsteps)
  for i=0, n_elements(pars[*,0])-1 do begin
     if tofit[3,i] eq -1 then begin
        pars[i,*] = str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value
        if str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).unit eq 'Radians' then angular = [angular,i]
     endif else begin
        pars[i,*] = (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].value
        if (*str.(tofit[0,i])[tofit[1,i]].(tofit[2,i])).(tofit[3,i])[tofit[4,i]].unit eq 'Radians' then angular = [angular,i]
     endelse
  endfor

  return, pars

end
