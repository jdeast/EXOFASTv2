function massradius_chenreverse, rp, tol=tol

  if n_elements(tol) eq 0 then tol = 1d-9

  minmp = 0
  maxmp = 1000

  while maxmp - minmp gt tol do begin

     mp = (maxmp+minmp)/2d0
     trialrp = massradius_chen(mp)
     
     if rp gt trialrp then minmp = mp $
     else maxmp = mp

;     print, minmp, maxmp, trialrp, rp

  endwhile



  return, mp
end
