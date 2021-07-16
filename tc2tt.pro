function tc2tt, time, e, inc, omega, period, tol=tol, reverse_correction=reverse_correction 

if n_elements(tol) eq 0 then tol = 1d-15

;; if inc=0, this will fail... but it should never transit
;bad = where(inc eq 0d0,nbad)
;if nbad eq 0 then return, time

;; if computing the reverse translation, the input time is T_T
if keyword_set(reverse_correction) then tt = time $
else tc = time

;; compute the phase at conjunction
phase=exofast_getphase(e,omega,/pri)

;; interatively compute theta, the angular offset between the 
;; time of conjunction and the time of minimum projected separation 
;; citation:
thetanew = 0d0
maxiter = 1d2
for i=0L, maxiter-1 do begin
   thetaold = thetanew
   thetanew = atan(-e*cos(omega)*cos(inc)^2/(cos(thetaold)*sin(inc)^2 + e*sin(omega)))
   if max(abs(thetanew - thetaold)) lt tol then break
endfor

if i ge maxiter then begin
   good = where(abs(thetanew - thetaold) lt tol,complement=bad)
   if bad[0] ne -1 then begin
      nbad = n_elements(bad)
      for j=0L, nbad-1 do begin
         print, time[bad[j]], e[bad[j]], inc[bad[j]], omega[bad[j]], $
                period[bad[j]], thetanew[bad[j]], thetaold[bad[j]],$
                (thetanew-thetaold)[bad[j]], format='(8(f0.20,x))'
      endfor
   endif
   message, "maximum iterations exceeded; please investigate"
endif

phase0 = exofast_getphase(e,omega,trueanom=!dpi/2d0 - omega + thetanew)

;; translate theta into an offset in time
if keyword_set(reverse_correction) then begin
   tc = tt - period*(phase0-phase)

   toohigh = where(tc gt tt+period/2d0)
   if toohigh[0] ne -1 then tc[toohigh] -= period
   toolow = where(tc lt tt-period/2d0)
   if toolow[0] ne -1 then tc[toolow] += period

   return, tc
endif else begin
   tt = tc + period*(phase0-phase)
   
   toohigh = where(tt gt tc+period/2d0)
   if toohigh[0] ne -1 then tt[toohigh] -= period
   toolow = where(tt lt tc-period/2d0)
   if toolow[0] ne -1 then tt[toolow] += period

   return, tt
endelse

end

 
