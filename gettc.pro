pro gettc, secondary=secondary

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

for i=0L, nplanets-1 do begin

   print
   print, 'planet ' + strtrim(i,2)

   for k=0, ntransits-1 do begin

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
