;; The Chen & Kipping relation is not monotonic!!! 
;; This function will never return masses between 
;; 0.278 and 0.414 M_Jupiter (88.49 to 131.58 M_Earth)
;; or between
;; 0.08 to 0.104 M_Sun (26644.8 to 34737.9 M_Earth, 83.83 to 109.3 M_Jupiter) 

function massradius_chenreverse, rp, tol=tol

  if rp lt 2.04d0^0.279d0 then begin
     mp = rp^(1d0/0.279d0)
  endif else if rp le 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.044d0)*26644.8321d0^(-0.044d0) then begin
     norm = 2.04d0^(0.279d0-0.589d0)
     mp = (rp/norm)^(1d0/0.589d0)
  endif else if rp le 2.04d0^(0.279d0-0.589d0)*131.58079d0^0.589d0 then begin
     norm = 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.044d0)
     mp = (rp/norm)^(-1d0/0.044d0)
  endif else begin
     norm = 2.04d0^(0.279d0-0.589d0)*131.58079d0^(0.589d0+0.044d0)*26644.8321d0^(-0.044d0-0.881d0)
     mp = (rp/norm)^(1d0/0.881d0)
  endelse

  return, mp

end
