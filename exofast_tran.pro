function exofast_tran, time, inc, ar, tp, period, e, omega, p, u1, u2, f0, $
                       rstar=rstar, thermal=thermal, reflect=reflect, $
                       dilute=dilute, tc=tc, q=q,x1=x1,y1=y1,z1=z1
AU = 215.094177d0

;; if we have the stellar radius, we can convert time to the
;; target's barycentric frame
if arg_present(rstar) then begin
   transitbjd = bjd2target(time, inclination=inc, a=ar*rstar, tp=tp, $
                           period=period, e=e, omega=omega,q=q)
endif else transitbjd = time

;; the impact parameter for each BJD
z = exofast_getb2(transitbjd, i=inc, a=ar, tperiastron=tp, period=period,$
                  e=e,omega=omega,z2=depth,x2=x,y2=y,q=q);,x0=x0,y0=y0,z0=z0)

;z = exofast_getb(transitbjd, i=inc, a=ar, tperiastron=tp, period=period,$
;                  e=e,omega=omega,z=depth,x=x,y=y)

;stop

if arg_present(z1) then depth += z1
if arg_present(x1) then x += x1
if arg_present(y1) then y += y1

;; Primary transit
modelflux = dblarr(n_elements(time))+1d0
primary = where(depth gt 0, complement=secondary)
if primary[0] ne - 1 then begin
   exofast_occultquad, z[primary], u1, u2, p, mu1
   modelflux[primary] =  mu1
endif

;; calculate the fraction of the planet that is visible for each time
if arg_present(thermal) or arg_present(reflect) then begin
   planetvisible = dblarr(n_elements(time)) + 1d0
   if secondary[0] ne - 1 then begin
      exofast_occultquad, z[secondary]/p, 0, 0, 1d0/p, mu1
      planetvisible[secondary] = mu1
   endif
endif

;; thermal emission from planet (isotropic)
if arg_present(thermal) then modelflux += 1d-6*thermal*planetvisible

;; phase-dependent reflection off planet
if arg_present(reflect) then $
   modelflux+=1d-6*reflect*cos(2d0*!dpi*(transitbjd-tc)/period)*planetvisible

;; dilution due to neighboring star
if arg_present(dilute) then modelflux *= ((1d0-dilute)+dilute)

;; normalization
modelflux *= f0

return, modelflux

end
