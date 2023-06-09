pro choplc, filename, tc, period, t14=t14, prefix=prefix, flatten=flatten

if n_elements(prefix) eq 0 then prefix = 'transit'

readcol, filename, time, flux, err, format='d,d,d'

if keyword_set(flatten) then begin
   norm = keplerspline(time, flux, breakp=breakp, ndays=0.75)
   flux /= norm
   err /= norm
endif


nplanets = n_elements(tc)

if n_elements(period) ne nplanets then message, 'Tc and Period must have the same number of elements'

if n_elements(t14) eq 0 then t14 = dblarr(nplanets) + 3d0/24d0 $
else if n_elements(t14) eq 1 then t14 = dblarr(nplanets) + t14 $
else if n_elements(t14) ne nplanets then message, 'Tc and Period must have the same number of elements'

planetletter = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

mintime = min(time,max=maxtime)


print, '***NOTE: Chopping the LC imposes these implicit uniform priors on the ephemeris.***'
print, 'The following priors should be added to the prior file to avoid artifacts from these implicit priors'
print

for i=0, n_elements(tc)-1 do begin

   minepoch = floor((mintime - tc[i])/period[i])
   maxepoch = ceil((maxtime - tc[i])/period[i])

   minepoch_real = !values.d_infinity
   maxepoch_real = -!values.d_infinity
   for j=minepoch, maxepoch do begin

      match = where(time ge (tc[i] + j*period[i] - t14[i]) and time le (tc[i] + j*period[i] + t14[i]))
      if match[0] ne -1 then begin
         if j lt minepoch_real then minepoch_real = j
         if j gt maxepoch_real then maxepoch_real = j

         caldat, tc[i]+j*period[i], month, day, year
         datestr = "n" + string(year,month,day,format='(i04,i02,i02)')

         forprint, time[match], flux[match], err[match], format='(f0.8,x,f0.8,x,f0.8)',/nocomment, textout=string(datestr,prefix,planetletter[i],j,'.dat', format='(a,a,a,".",i04,a)'),/silent
      endif

   endfor

   if finite(minepoch_real) and finite(maxepoch_real) then begin
      ;; chopping the LC like this imposes these implicit priors on the model
      period_error = t14[i]/(maxepoch_real-minepoch_real)
      print, string(i,tc[i],tc[i]-t14[i],tc[i]+t14[i],format='("tc_",i1,x,f0.6," -1 ",f0.6,x,f0.6)')
      print, string(i,period[i],period[i]-period_error,period[i]+period_error,format='("period_",i1,x,f0.12," -1 ",f0.12,x,f0.12)')
   endif

endfor

end
