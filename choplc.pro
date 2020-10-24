pro choplc, filename, tc, period, t14=t14, prefix=prefix

if n_elements(prefix) eq 0 then prefix = 'transit'

readcol, filename, time, flux, err, format='d,d,d'

nplanets = n_elements(tc)

if n_elements(period) ne nplanets then message, 'Tc and Period must have the same number of elements'

if n_elements(t14) eq 0 then t14 = dblarr(nplanets) + 3d0/24d0 $
else if n_elements(t14) eq 1 then t14 = dblarr(nplanets) + t14 $
else if n_elements(t14) ne nplanets then message, 'Tc and Period must have the same number of elements'

planetletter = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

mintime = min(time,max=maxtime)

for i=0, n_elements(tc)-1 do begin

   minepoch = floor((mintime - tc)/period)
   maxepoch = ceil((maxtime - tc)/period)
   nepochs = maxepoch-minepoch+1
   epochs = minepoch + lindgen(nepochs)

   for j=0, nepochs-1 do begin

      match = where(time ge (tc[i] + j*period[i] - t14[i]) and time le (tc[i] + j*period[i] + t14[i]))
      if match[0] ne -1 then begin
         forprint, time[match], flux[match], err[match], format='(f0.8,x,f0.8,x,f0.8)',/nocomment, textout=string(prefix,planetletter[i],j,'.dat', format='(a,a,i04,a)')
      endif
   endfor
endfor

end
