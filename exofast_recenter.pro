;+
; NAME:
;   EXOFAST_RECENTER
;
; PURPOSE: 
;   Recenters a distribution of periodic parameters to the domain
;      (mode-period/2,mode+period/2]
;
; CALLING SEQUENCE:
;   par = recenter(par,period)
;
; INPUTS:
;   PAR    = An array of parameters
;   PERIOD = The period of the parameter distribution
;
; EXAMPLE:
;
; MODIFICATION HISTORY
;  2012/06 -- Jason Eastman (LCOGT)
;-

function exofast_recenter, par, period, logname=logname

hist = histogram(par,nbins=100,locations=x)
max = max(hist,modendx)
mode = x[modendx]

if n_elements(period) eq 1 then per = replicate(period,n_elements(par)) $
else if n_elements(period) eq n_elements(par) then per = period $
else begin
   printandlog, "period must have 1 or npar elements", logname
   stop
endelse


bad = where(mode - period/2d0 - mode eq 0d0, nbad)
if nbad ne 0 then begin
   printandlog, 'ERROR: the period is negligible compared to the value', logname
   printandlog, 'This usually means the transit time ran away and requires some bound', logname
   printandlog, 'Consider a wide, uniform prior on Tc and/or set REJECTFLATMODEL', logname
   printandlog, 'Returning without recentering the distribution', logname
   printandlog, '***THIS RESULT SHOULD ONLY BE USED TO DEBUG THE PROBLEM!***', logname
   return, par
endif

nper = round((mode - par)/per,/L64)
par -= nper*per

repeat begin
   toohigh = where(par gt (mode + period/2d0))
   if toohigh[0] ne -1 then par[toohigh] -= per[toohigh]
endrep until toohigh[0] eq -1

repeat begin
   toolow = where(par le (mode - period/2d0))
   if toolow[0] ne -1 then par[toolow] += period[toolow]
endrep until toolow[0] eq -1

return, par

end
