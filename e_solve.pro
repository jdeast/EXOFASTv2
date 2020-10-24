;; compute T_FWHM (~mid ingress to mid egress)
function tfwhm, period, ar, inc, e, omega

  cosi = cos(inc)
  sini = sin(inc)
  esinw = e*cos(omega)
  
  b = ar*cosi*sqrt(1d0-e^2)/(1d0+esinw)
  return, period/!dpi*asin(sqrt(1d0 - b^2)/(sini*ar))*sqrt(1d0-e^2)/(1d0+esinw)
  
end
  
function e_root, e

  common e_block, tfwhm0, ar, inc, omega
  
  tfwhmc = tfwhm(1d0, ar, inc, 0d0, omega)
  tfwhme = tfwhm(1d0, ar, inc,   e, omega)

  return, tfwhme/tfwhmc - tfwhm0
  
end

;; there are 0-2 roots, which may be complex
function e_solve, tfwhm_0, ar_0, inc_0, omega_0, solution

  tol = 1d-9

  common e_block, tfwhm0, ar, inc, omega
  tfwhm0 = tfwhm_0
  ar = ar_0
  inc = inc_0
  omega = omega_0

if 0 then begin
  ;period = 5d0
  ;ar = 11d0
  ;inc = 89.7d0*!dpi/180d0
  ;omega = -90d0*!dpi/180d0
  ;tfwhm0 = 2.1d0

  mine = -1d0
  maxe = 1d0
  nsteps = 1d4
  e = mine + (maxe-mine)*dindgen(nsteps)/(nsteps-1)
  tfwhm0 = tfwhm(period, ar, inc, e, omega)
  plot, e, tfwhm0
;stop
endif


  eguess1c = [complex(0d0,-0.5d0), complex(0d0,-0.75d0), complex(0d0,-1d0)]
  eguess2c = [complex(0d0,0d0)   , complex(0d0,-0.25d0), complex(0d0,-0.5d0)]
  eguess3c = [complex(0d0,0d0)   , complex(0d0,0.25d0), complex(0d0,0.5d0)]
  eguess4c = [complex(0d0,0.5d0) , complex(0d0,0.75d0), complex(0d0,1d0)]

  eguess1 = [complex(-0.5d0,0d0), complex(-0.75d0,0d0), complex(-1d0,0d0)]
  eguess2 = [complex(0d0,0d0)   , complex(-0.25d0,0d0), complex(-0.5d0,0d0)]
  eguess3 = [complex(0d0,0d0)   , complex(0.25d0,0d0) , complex(0.5d0,0d0)]
  eguess4 = [complex(0.5d0,0d0) , complex(0.75d0,0d0) , complex(1d0,0d0)]

;  t0 = systime(/seconds)
  e1 = fx_root(eguess1, 'e_root', /double, tol=tol)
  e2 = fx_root(eguess2, 'e_root', /double, tol=tol)
  e3 = fx_root(eguess3, 'e_root', /double, tol=tol)
  e4 = fx_root(eguess4, 'e_root', /double, tol=tol)
  e1c = fx_root(eguess1c, 'e_root', /double, tol=tol)
  e2c = fx_root(eguess2c, 'e_root', /double, tol=tol)
  e3c = fx_root(eguess3c, 'e_root', /double, tol=tol)
  e4c = fx_root(eguess4c, 'e_root', /double, tol=tol)
;  print, systime(/seconds)-t0

  e = round(real_part([e1,e2,e3,e4,e1c,e2c,e3c,e4c])*1d9,/l64)/1d9
  e = e[uniq(e)]

  bad = where(abs(tfwhm(1d0, ar, inc, e, omega)/tfwhm(1d0, ar, inc, 0, omega)- tfwhm0) lt tol)
  if bad[0] ne -1 then e[bad] = !values.d_nan

  print, e

  if n_elements(e) gt 1 and solution then return, e[1]
  return, e[0]

  print
  print, e
  print
  print, tfwhm(period, ar, inc, e, omega)/tfwhm(period, ar, inc, 0, omega), tfwhm0

;
;  print, real_part(e1), real_part(e_root(e1)), real_part(e2), real_part(e_root(e2))

stop

  print, e

end

function test

bc = cosi/ar
tfwhmc = period/!dpi*asin(sqrt(1d0 - bc^2)/(sin(inc)*ar))

nsteps = 1d3
e = dindgen(nsteps)/(nsteps-1)

emin = 0d0
emax = 1d0

tol = 1d-9
repeat begin

   e = (emin+emax)/2d0
   b = cosi/ar*sqrt(1d0-e^2)/(1d0+e*sin(omega))
   tfwhm = period/!dpi*asin(sqrt(1d0 - b^2)/(sin(inc)*ar))*sqrt(1d0-e^2)/(1d0+e*sin(omega))/tfwhmc

   
  
   if tfwhm gt tfwhm0 then emin = e $
   else emax = e

   print, e, emin, emax, tfwhm, tfwhm0

endrep until emax-emin lt tol


stop
if abs(tfwhm - tfwhm0) gt 1d-6 then return, !values.nan
return, e


;plot, e, tfwhm
;oplot, [-9d9,9d9],[1d0/24d0,1d0/24d0]


stop
end
