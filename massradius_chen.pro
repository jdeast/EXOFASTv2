function massradius_chen, mp, rperr=rperr

  if mp lt 2.04d0 then begin
     rp = mp^0.279d0
     rperr = rp*0.0403d0
  endif else if mp lt 131.58079d0  then begin ;; 0.414 M_jupiter
     norm = 2.04d0^(0.279d0-0.589d0)
     rp = norm*mp^0.589d0
     rperr = rp*0.1460d0
  endif else begin
     norm = 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.44d0)
     rp = norm*mp^(-0.44d0)
     rperr = rp*0.0737d0
  endelse

  return, rp

end
