function massradius_chen, mp, rperr=rperr

  if mp le 2.04d0 then begin
     rp = mp^0.279d0
     rperr = rp*0.0403d0
  endif else if mp le 131.58079d0  then begin ;; 0.414 M_jupiter
     norm = 2.04d0^(0.279d0-0.589d0)
     rp = norm*mp^0.589d0
     rperr = rp*0.1460d0
  endif else if mp le 26644.8321d0 then begin ;; 0.08 M_sun
     norm = 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.044d0)
     rp = norm*mp^(-0.044d0)
     rperr = rp*0.0737d0
  endif else begin
     norm = 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.044d0)*26644.8321d0^(-0.044d0-0.881d0)
     rp = norm*mp^(0.881d0)
     rperr = 0.043d0     
  endelse

  return, rp

end
