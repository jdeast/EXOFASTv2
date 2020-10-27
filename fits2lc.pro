pro fits2lc, fitsname, t14=t14, tc=tc, period=period, flatten=flatten, splinespace=splinespace

if n_elements(splinespace) eq 0 then splinespace = 0.75d0

nplanets = n_elements(tc)

fileparts = strsplit(file_basename(fitsname,'.fits'),'-',/extract)
sector = string(long(strsplit(fileparts[1],'s',/extract)),format='(i02)')
ticid = fileparts[2]

;; from SPOC (2 minute cadence)
data = mrdfits(fitsname,1)
time = data.time + 2457000d0 ;; BJD_TDB?

flux = data.pdcsap_flux
fluxerr = data.pdcsap_flux_err

;; find the momentum dumps to introduce a break in the spline fit
;; (where quality bit 6 is set)
momentumdump = bytarr(n_elements(data.quality))
for k=0L, n_elements(data.quality)-1 do $
   momentumdump[k] = (bitmask2arr((data.quality)[k],nbits=12))[5]
breaktimes = time[where(momentumdump)]

;; remove flagged bad data points
good = where(finite(time) and time ne 0 and finite(flux) and data.quality eq 0L)
time = time[good]
flux = flux[good]
fluxerr = fluxerr[good]

;; now define break points for each momentum dump
breakp = [-1]
for k=0L, n_elements(breaktimes)-1 do begin
   arr = time-breaktimes[k]
   neg = where(arr lt 0)
   if neg[0] ne -1 then arr[neg] = !values.d_infinity
   junk = min(arr,ndx)
   breakp = [breakp,ndx]
endfor

;; remove repeats
sorted = sort(breakp)
breakp = breakp[sorted]
breakp = breakp[uniq(breakp)]
;; remove -1 if more than 1 element
if n_elements(breakp) gt 1 then breakp = breakp[1:n_elements(breakp)-1]

;; if Tc, period, and T14 are not given, use all points
if n_elements(period) eq 0 or n_elements(tc) eq 0 or n_elements(t14) eq 0 then begin
   use = lindgen(n_elements(time))
   intransit = [1]
endif else begin
   ;; otherwise, trim to +/- 2*T14
   phase = dblarr(nplanets,n_elements(time))
   use = []
   intransit = []
   for k=0L, nplanets-1 do begin
      phase[k,*] = (time-tc[k]) mod period[k]
      toohigh = where(phase[k,*] gt period[k]/2d0)
      if toohigh[0] ne -1 then phase[k,toohigh] -= period[k]
      toolow = where(phase[k,*] le -period[k]/2d0)
      if toolow[0] ne -1 then phase[k,toolow] += period[k]
      thisuse = where(abs(phase[k,*]) le (2d0*duration[k]))
      if thisuse[0] ne -1 then use = [use,thisuse]
      thisintransit = where(phase[k,*] gt (-duration[k]/2d0) and phase[k,*] lt (duration[k]/2))
      if thisintransit[0] ne -1 then intransit = [intransit,thisintransit]
   endfor
   
   ;; remove overlapping elements
   if n_elements(use) gt 0 then begin
      sorted = sort(use)
      use = use[sorted]
      use = use[uniq(use)]
   endif
endelse

if n_elements(use) gt 3 and n_elements(intransit) ge 1 then begin
   ;; flatten light curve
   if keyword_set(flatten) then begin
      norm = keplerspline(time, flux, breakp=breakp, ndays=splinespace)
   endif else begin
      norm = mean(flux)
   endelse

   caldat, time[use[0]], month, day, year
   datestr = string(year,month,day,format='("n",i04,i02,i02)')
   
   tranfile = datestr + '.TESS.TESS.s' + sector + '.spoc.dat'
   
   exofast_forprint, time[use],(flux/norm)[use],(fluxerr/norm)[use],format='(f0.6,x,f0.6,x,f0.6)',textout=tranfile,/nocomment,/silent
endif else begin
   print, tranfile + ' has no transit; skipping'
endelse

end
