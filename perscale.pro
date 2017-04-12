pro perscale, ss

for i=0, ss.nplanets-1 do begin

   ntransits = 0L
   allminepoch = !values.d_infinity
   allmaxepoch = -!values.d_infinity
   for j=0, ss.ntran-1 do begin
      
      transit = *(ss.transit[j].transitptrs)
      mintime = min(transit.bjd,max=maxtime)
      
      epoch = (transit.bjd - ss.planet[i].tc.value)/ss.planet[i].period.value
      minepoch = fix(min(epoch))
      maxepoch = fix(max(epoch))

      ntransits += maxepoch - minepoch + 1
      if minepoch lt allminepoch then allminepoch = minepoch
      if maxepoch gt allmaxepoch then allmaxepoch = maxepoch

   endfor
   
   ;; we can do better than this (Yee & Gaudi, 2008; Cummings 2006)
   if (allmaxepoch-allminepoch+1) gt 1 then begin
      ;; with multiple transits, the period is well known
      uper = sqrt(2d0)*0.01/sqrt(allmaxepoch-allminepoch)
      ulogp = uper/(alog(10d0)*ss.planet[i].period.value)
      ss.planet[i].logp.scale = 3d0*ulogp

;      ss.planet[i].logp.scale = 0.001
   endif else if ss.planet[i].fitrv then begin
      ;; with only RVs, the period is not as well known
      uper = ss.planet[i].period.value/100d0
      ulogp = uper/(alog(10d0)*ss.planet[i].period.value)
      ss.planet[i].logp.scale = ulogp
   endif else begin
      ;; if we have a single transit and no RVs, the period degenerate
      ;; with eccentricity and highly uncertain
      ss.planet[i].logp.scale = 1.5
   endelse

endfor


end
