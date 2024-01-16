function tc2tt, time, e, inc, omega, period, tol=tol, tt2tc=tt2tc, ts2te=ts2te, te2ts=te2ts

if n_elements(tol) eq 0 then tol = 1d-15

;; if inc=0, this will fail... but it should never transit
;bad = where(inc eq 0d0,nbad)
;if nbad eq 0 then return, time

;; if computing the reverse translation, the input time is T_T
;if keyword_set(tt2tc) then tt = time $
;else if keyword_set(ts2te) then ts = time $
;else if keyword_set(te2ts) then te = time $
;else tc = time

;; interatively compute theta, the angular offset between the 
;; time of conjunction and the time of minimum projected separation 
;; see equation 14 of Csizmadia, 2020; Marynov 1973;
thetanew = 0d0
maxiter = 1d2
for i=0L, maxiter-1 do begin
   thetaold = thetanew

   if keyword_set(ts2te) or keyword_set(te2ts) then begin
      ;; computing secondary eclipse time
      thetanew = -atan(-e*cos(omega)*cos(inc)^2/(cos(thetaold)*sin(inc)^2 - e*sin(omega)))
   endif else begin  
      ;; computing primary transit time
      thetanew = atan(-e*cos(omega)*cos(inc)^2/(cos(thetaold)*sin(inc)^2 + e*sin(omega)))
   endelse

   if max(abs(thetanew - thetaold)) lt tol then break
endfor

if i ge maxiter then begin
   good = where(abs(thetanew - thetaold) lt tol,complement=bad)
   if bad[0] ne -1 then begin
      nbad = n_elements(bad)
      print, "these parameter combinations are not converging after " + strtrim(maxiter,2) + " iterations:"
      
      print, '       BJD_TDB                          e                 inc (radians)          omega (radians)         period (days)          theta_new              theta_old                diff'
      for j=0L, nbad-1 do begin
         print, time[bad[j]], e[bad[j]], inc[bad[j]], omega[bad[j]], $
                period[bad[j]], thetanew[bad[j]], thetaold[bad[j]],$
                (thetanew-thetaold)[bad[j]], format='(8(f0.20,x))'
      endfor
   endif
   message, "maximum iterations exceeded; please investigate"
endif

if keyword_set(ts2te) or keyword_set(te2ts) then begin
   ;; compute the phase at superior conjunction
   phase=exofast_getphase(e,omega,/sec)
   phase0 = exofast_getphase(e,omega,trueanom=3d0*!dpi/2d0 - omega + thetanew)
   if keyword_set(te2ts) then begin 
      ;; computing te from ts
;      te = ts - period*(phase0-phase)
      tfinal = time - period*(phase0-phase)
   endif else begin
      ;; computing ts from te
;      ts = te + period*(phase0-phase)
      tfinal = time + period*(phase0-phase)
   endelse
endif else begin
   ;; compute the phase at inferior conjunction
   phase=exofast_getphase(e,omega,/pri)
   phase0 = exofast_getphase(e,omega,trueanom=!dpi/2d0 - omega + thetanew)
   if keyword_set(tt2tc) then begin
      ;; computing tt from tc
;      tc = tt - period*(phase0-phase)
      tfinal = tt - period*(phase0-phase)
   endif else begin
      ;; computing tc from tt
;      tt = tc + period*(phase0-phase)
      tfinal = time + period*(phase0-phase)
   endelse
endelse

;; make sure the output is within +/- period of the input time
toohigh = where(tfinal gt time+period/2d0)
if toohigh[0] ne -1 then tfinal[toohigh] -= period
toolow = where(tfinal lt time-period/2d0)
if toolow[0] ne -1 then tfinal[toolow] += period

return, tfinal

end

 
