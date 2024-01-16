; COMMON BLOCKS:
;   RV_BLOCK - See exofast.pro for definition
;
; MODIFICATION HISTORY 
;  2023/09 -- Jason Eastman (CfA)
;-
function exofast_getchi2_rv_fitcir, pars,time,rv,err,slope=slope,quad=quad

; pars[0] = time of transit center
; pars[1] = period
; pars[2] = e*cos(omega)
; pars[3] = e*sin(omega)
; pars[4] = K
; pars[5] = gamma
; pars[6] = slope
; pars[7] = quad

;; as defined by circular orbit
pars[2] = 0.d0
pars[3] = 0.d0
e = 0
omega = !dpi/2.d0

mintime = min(time,max=maxtime)
t0 = (mintime + maxtime)/2.d0

;; fit the amplitude, phase, offset, and slope analytically 
derivs = transpose([[cos(2.d0*!dpi*time/pars[1])/err],$
                    [sin(2.d0*!dpi*time/pars[1])/err],$
                    [1.d0/err]])

if keyword_set(slope) or keyword_set(quad) then $
   derivs = [derivs,transpose((time-t0)/err)]

if keyword_set(quad) then $
   derivs = [derivs,transpose((time-t0)^2/err)]

;; the magic
npars = n_elements(derivs[*,0])
datarr = replicate(1,npars)#rv
errarr = replicate(1,npars)#err
b = matrix_multiply(derivs,derivs,/btranspose)
d = total(derivs*datarr/errarr,2)
a = invert(b)#d

K = sqrt(a[0]^2 + a[1]^2)
phi = -atan(a[1]/a[0])
if a[0] lt 0 then phi += !dpi

pars[0] = pars[1]*(omega-phi)/(2.d0*!dpi)
pars[4] = K
pars[5] = a[2]

;; calculate the model
model = K*cos(2*!dpi*time/pars[1] + phi)+pars[5]

if keyword_set(slope) or keyword_set(quad) then begin
   pars[6] = a[3]
   model += pars[6]*(time-t0)
endif
if keyword_set(quad) then begin
   pars[7] = a[4]
   model += pars[7]*(time-t0)^2
endif

;; compute the chi2
chi2 = total(((rv - model)/err)^2)

return, chi2

end

;+
; NAME:
;   EXOFAST_LOMBSCARGLE
;
; PURPOSE: 
;   This function returns the NMIN best periods determined by a
;   lomb-scargle periodogram. 
; 
; CALLING SEQUENCE:
;   periods = exofast_lombscargle(time, rv [BESTPARS=,
;   SCALE=, NMIN=, MINPERIOD=, MAXPERIOD=,NOSLOPE=, /NYQUIST, /PLOT])
;
; INPUTS:
;   TIME        - A time (eg, BJD_TDB) for each of the RV data points
;   RV          - The RV data points
;
; OPTIONAL OUTPUTS:  
;   BESTPARS    - The best (analytically fit) parameters determined from
;                 the lomb-scargle periodogram (assumes circular orbit).
;   SCALE       - The scale (~error) of the period, determined by the
;                 spacing of the sampling of the periodogram.
;   NMIN        - The number of periods corresponding to chi^2 minima
;                 to return (default=5)
;   MINPERIOD   - The minimum period to explore. Default = 1 unit of
;                 time (typically days). This value is not used if
;                 /NYQUIST is set.
;   MAXPERIOD   - The maximum period to explore. Default is the full
;                 range of dates.
;
; OPTIONAL KEYWORDS
;   SLOPE       - If set, a slope trend is analytically fit to the RV data.
;   QUAD        - If set, a quadratic trend is analytically fit to the RV data.
;   NYQUIST     - If set, the minimum period explored is the nyquist
;                 frequency (half the minimum spacing between two
;                 times). This takes precedence over MINPERIOD. Unless
;                 the times are uniformly spaced (which is not
;                 typical), use of this keyword is not recommended.
;   PLOT        - A diagnostic tool. If set, the lomb-scargle
;                 periodogram will be plotted, with the best period(s)
;                 overplotted. This is useful for determining good
;                 bounds, or if the observed behavior is unexpected.
;
; DEPENDENCIES:
;   BUIELIB (http://www.boulder.swri.edu/~buie/idl/)
;
; OUTPUTS:
;    result     - An NMIN elements array containing the best periods
;                 determined from the lomb-scargle periodogram.
;
; MODIFICATION HISTORY 
;  2023/09 -- Jason Eastman (Cfa)
;-

function exofast_lombscargle,time,rv,err,bestpars=bestpars,scale=scale,nmin=nmin,$
                             minperiod=minperiod,maxperiod=maxperiod,$
                             slope=slope,quad=quad,plot=plot, periods=periods, chi2=chi2,psname=psname
                              

mintime = min(time,max=maxtime)

if n_elements(nmin) eq 0 then nmin = 5
if n_elements(minperiod) eq 0 then minperiod = 1.d0
if n_elements(maxperiod) eq 0 then maxperiod = (maxtime - mintime)

;; optimal period sampling
duration = max(time) - min(time)
periods = double(minperiod)
np = 1L
repeat begin
   periods = [periods,periods[np-1L] + periods[np-1L]^2/(4d0*!dpi*duration)]
   np++
endrep until periods[np-1L] ge maxperiod

chi2 = dblarr(np)
pars = dblarr(7,np)

minchi2 = !values.d_infinity

for i=0L, np-1 do begin
    initpars = [0d0,periods[i],0d0,0d0,0d0,0d0,0d0,0d0]
    chi2[i] = exofast_getchi2_rv_fitcir(initpars,time,rv,err,slope=slope,quad=quad)
    pars[*,i] = initpars
endfor

;; find the NMIN local minima
minima = lclxtrem(chi2)
nminima = n_elements(minima)
best = minima[nminima - nmin > 0:nminima-1]

bestpars = pars[*,best]
scale = periods[best+1] - periods[best-1]

;; plot the lomb-scargle periodogram
if keyword_set(plot) then begin
   if keyword_set(psname) then begin
      mydevice = !d.name
      set_plot, 'PS'
      device, filename=psname,/color,bits=24
      loadct, 39, /silent
      red = 254
   endif else begin
      window, 0, retain=2
      red = '0000ff'x
   endelse
   
    plot, periods, chi2, xtitle='Period (days)',ytitle=textoidl('\chi^2'),/xlog,/xstyle,/ystyle
    oplot, periods[best],chi2[best],color=red,psym=1

    if keyword_set(psname) then begin
       device, /close
       set_plot, mydevice
    endif else begin
       print, 'Type ".con" to continue'
       stop
    endelse
endif

return, periods[best]

end
