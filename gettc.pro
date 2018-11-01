;+
;
; MCMCSS - The MCMC Stellar structure output by exofastv2.pro 
;
; EPOCHS - An NPLANETS pointer array (must match mcmcss.nplanets),
;          pointing to an array specifying the epochs of the desired
;          transit times. A null pointer signifies no output desired
;          for that planet.
;
; NSAMPLES - This calculation is expensive and may take a decent
;            fraction of the total MCMC runtime. We randomly sample
;            NSAMPLES from the MCMCSS values to compute the
;            uncertainty in each transit time. Default is 1000. Use -1
;            or a value larger than mcmcss.nvalues to compute the Tcs
;            at each epoch for all values.
;
; RESULT - An NPLANETS pointer array, pointing to an NEPOCHS x NSAMPLES
;          array of Transit times, including stellar reflex motion and
;          fitted TTVs.
; 
;-
function gettc, mcmcss, epochs, nsamples=nsamples, secondary=secondary, tol=tol

if n_elements(tol) eq 0 then tol = 0.1d0/86400d0 ;; accurate to 0.1 seconds

chi2 = reform((*(mcmcss.chi2)),mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains)
ngoodlinks = n_elements(goodchains)*(mcmcss.nsteps/mcmcss.nchains - burnndx)

if n_elements(nsamples) eq 0 then nsamples = 10000d0
if nsamples eq -1 then nsamples = ngoodlinks

;; choose NSAMPLES random links from random chains
if nsamples ge ngoodlinks then ndx = lindgen(ngoodlinks) $
else ndx = randomu(seed,nsamples)*ngoodlinks

;tcs = ptrarr(mcmcss.nplanets,/allocate_heap)
;for j=0L, mcmcss.nplanets-1L do *(tcs[j]) = dblarr(mcmcss.ntran,nsamples)

tcs = dblarr(mcmcss.nplanets,mcmcss.ntran,nsamples)
epochs = lonarr(mcmcss.nplanets,mcmcss.ntran)

for j=0L, mcmcss.nplanets-1 do begin
   
   allperiods = (reform((reform(mcmcss.planet[j].period.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains],ngoodlinks))
   alltcs = (reform((reform(mcmcss.planet[j].tc.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains],ngoodlinks))
   allt14s = (reform((reform(mcmcss.planet[j].t14.value,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains))[burnndx:*,goodchains],ngoodlinks))
   

   for i=0L, nsamples-1 do begin

   ;; get the motion of the star due to all planets
;   junk = exofast_getb2(transitbjd,inc=mcmcss.planet.i.value[ndx[i]],a=mcmcss.planet.ar.value[ndx[i]],$
;                        tperiastron=mcmcss.planet.tp.value[ndx[i]],$
;                        period=mcmcss.planet.period.value[ndx[i]],$
;                        e=mcmcss.planet.e.value[ndx[i]],omega=mcmcss.planet.omega.value[ndx[i]],$
;                        q=mcmcss.star.mstar.value[ndx[i]]/mcmcss.planet.mpsun.value[ndx[i]],$
;                        x1=x1,y1=y1,z1=z1)

      period = allperiods[ndx[i]]
      tc = alltcs[ndx[i]]
      t14 = allt14s[ndx[i]]
      
      for k=0L, mcmcss.ntran-1 do begin
         
         minbjd = min((*(mcmcss.transit[k].transitptrs)).bjd,max=maxbjd)
         epoch = round(((maxbjd+minbjd)/2d0 - tc)/period)

         thistc = tc + epoch*period
         if (minbjd le thistc + t14 and maxbjd ge thistc - t14) then begin
            
            
            maxbjd = tc + period*epoch + t14
            minbjd = tc + period*epoch - t14
            t0 = floor(minbjd)

            repeat begin
               
               transitbjd = (maxbjd + minbjd)/2d0

               z = exofast_getb2(transitbjd,inc=mcmcss.planet.i.value[ndx[i]],a=mcmcss.planet.ar.value[ndx[i]],$
                                 tperiastron=mcmcss.planet.tp.value[ndx[i]],$
                                 period=mcmcss.planet.period.value[ndx[i]],$
                                 e=mcmcss.planet.e.value[ndx[i]],omega=mcmcss.planet.omega.value[ndx[i]],$
                                 q=mcmcss.star.mstar.value[ndx[i]]/mcmcss.planet.mpsun.value[ndx[i]],$
                                 x0=x0,y0=y0,z0=z0)

               if x0[j] gt 0 then maxbjd = transitbjd $
               else minbjd = transitbjd

;               print, x0[j], minbjd-t0, maxbjd-t0, (maxbjd-minbjd)*86400d0, i, j, k;minbjd-t0, maxbjd-t0, transitbjd-t0

;               wait, 1
;               stop

            endrep until abs(maxbjd-minbjd) lt tol

            tcs[j,k,i] = transitbjd + mcmcss.transit[k].ttv.value[ndx[i]]
            epochs[j,k] = epoch
         endif

      endfor
   endfor
endfor   

planetname = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

for j=0L, mcmcss.nplanets-1 do begin
   for k=0L, mcmcss.ntran-1 do begin
      bad = where(tcs[j,k,*] eq 0d0)
      if bad[0] eq -1 then begin
         print, planetname[j], epochs[j,k], mean(tcs[j,k,*]), stddev(tcs[j,k,*]),mcmcss.transit[k].label, format='(a1,x,i3,x,f0.5,x,f0.5,x,a)'
      endif
   endfor
endfor


stop

return, tcs

stop


npoints = 1d5

G = 2942.71377d0 ;; R_sun^3/(m_sun*day^2), Torres 2010
nplanets = 3d0
period = [1d0,2d0,10d0]
period = 1+randomu(seed,3)*30d0
period = period[sort(period)]
tc = [0d0,3d0,5d0]
tc = randomu(seed,3)*max(period)
inc=[!dpi/2d0,!dpi/2d0,!dpi/2d0]
m1 = 1d0
m2 = [0.001d0,0.01d0,0.01d0]
q = m1/m2
a = (period^2*G*(m1+m2)/(4d0*!dpi^2))^(1d0/3d0)
;e = [0.5d0,0.25d0,0.9d0]
e = dblarr(3);[0.5d0,0.25d0,0.9d0]
omega = [!dpi/5d0,!dpi/6d0,!dpi/7d0]


bjd = min(tc) + dindgen(npoints)/(npoints-1)*max(period)*10
tp = dblarr(nplanets)
for i=0, nplanets-1L do begin
   phase = exofast_getphase(e[i],omega[i], /primary)
   tp[i] = tc[i] - period[i]*phase
endfor

if 0 then begin
b = exofast_getb2(bjd, inc=inc,a=a,tperiastron=tp, period=period,e=e,omega=omega,q=q,x0=x,y0=y,z0=z,x1=x1,y1=y1,z1=z1,x2=x2,y2=y2,z2=z2)

for i=0, nplanets-1 do begin
      
   ndx = lclxtrem(b[i,*])
!p.multi=[0,2,2]
   plot, b[i,*]
   oplot, b[i,ndx], b[i,ndx], psym=1, color='0000ff'x

endfor
endif

a2 = a*q/(1d0+q)
a1 = a2/q  
tol = 1d-8 ;; 1 ms
ntransits = 10



;; get the motion of the star due to all planets
junk = exofast_getb2(transitbjd,inc=ss.planet.i.value,a=ss.planet.ar.value,$
                     tperiastron=ss.planet.tp.value,$
                     period=ss.planet.period.value,$
                     e=ss.planet.e.value,omega=ss.planet.omega.value,$
                     q=ss.star.mstar.value/ss.planet.mpsun.value,$
                     x1=x1,y1=y1,z1=z1)


for i=0L, mcmcss.nplanets-1 do begin

   for k=0, mcmcss.ntran-1 do begin

      t0 = systime(/seconds) 
      minbjd = tc[i]-period[i]/4d0+k*period[i]
      maxbjd = tc[i]+period[i]/4d0+k*period[i]
      repeat begin
         
         bjd = (minbjd + maxbjd)/2d0
         
         x1 = 0d0
         z1 = 0d0
         for j=0L, nplanets-1 do begin
            meananom = (2.d0*!dpi*(1.d0 + (bjd - tp[j])/period[j])) mod (2.d0*!dpi)
            eccanom = exofast_keplereq(meananom, e[j])
            trueanom = 2d0*atan(sqrt((1d0 + e[j])/(1d0 - e[j]))*tan(0.5d0*eccanom))
            
            ;; calculate the star position in the barycentric frame
            r1 = a1[j]*(1d0-e[j]^2)/(1d0+e[j]*cos(trueanom))
            
            ;; rotate to observer's plane of reference
            x1 += -r1*cos(trueanom + omega[j] + !dpi)
            z1 += r1*sin(trueanom + omega[j] + !dpi)*sin(inc[j])
         endfor
         
         meananom = (2.d0*!dpi*(1.d0 + (bjd - tp[i])/period[i])) mod (2.d0*!dpi)
         eccanom = exofast_keplereq(meananom, e[i])
         trueanom = 2d0*atan(sqrt((1d0 + e[i])/(1d0 - e[i]))*tan(0.5d0*eccanom))
         
         r = a2[i]*(1d0-e[i]^2)/(1d0+e[i]*cos(trueanom))
         x = -r*cos(trueanom + omega[i])-x1
         z = r*sin(trueanom + omega[i] + !dpi)*sin(inc[i]) -z1

         ;; bracket the primary or secondary
         if keyword_set(secondary) then begin
            if z le 0d0 and x gt 0d0 then minbjd = bjd $
            else if z le 0d0 and x le 0d0  then maxbjd = bjd $
            else if z gt 0d0 and x gt 0d0 then maxbjd = bjd $
            else if z gt 0d0 and x le 0d0 then minbjd = bjd 
         endif else begin
            if z le 0d0 and x gt 0d0 then maxbjd = bjd $
            else if z le 0d0 and x le 0d0  then minbjd = bjd $
            else if z gt 0d0 and x gt 0d0 then minbjd = bjd $
            else if z gt 0d0 and x le 0d0 then maxbjd = bjd 
         endelse
         
         if x gt 0 then begin
            maxbjd = bjd
         endif else minbjd = bjd
         
      endrep until abs(x) le tol
      
      print, bjd, (bjd-tc[i]-k*period[i])*86400d0, systime(/seconds)-t0
   endfor
endfor


end
