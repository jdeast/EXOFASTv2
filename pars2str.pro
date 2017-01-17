pro pars2str, pars, str

  tofit = *(str.tofit)
  npars = n_elements(tofit[0,*])
  if (size(pars))[0] eq 1 then begin
     ;; populate the structure with the parameter array
     for i=0, npars-1 do $
        str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = pars[i]
  endif else begin
     for i=0, npars-1 do $
        str.(tofit[0,i])[tofit[1,i]].(tofit[2,i]).value = transpose(pars[i,*])
  endelse     

end
