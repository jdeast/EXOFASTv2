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

function exofast_recenter, par, period

hist = histogram(par,nbins=100,locations=x)
max = max(hist,modendx)
mode = x[modendx]

if n_elements(period) eq 1 then per = replicate(period,n_elements(par)) $
else if n_elements(period) eq n_elements(par) then per = period $
else message, "period must have 1 or npar elements"

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
