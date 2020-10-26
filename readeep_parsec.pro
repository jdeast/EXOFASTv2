function readeep_parsec, mstar, initfeh, vvcrit=vvcrit, alpha=alpha, rstar=rstar, teff=teff, age=age,zxsun=zxsun

if n_elements(zxsun) eq 0 then zxsun=0.0207d0

if n_elements(alpha) eq 0 then alpha = 0d0
if n_elements(vvcrit) eq 0 then vvcrit = 0d0

;; check for allowed [Fe/H]
;; z grid points (these must be in order of [Fe/H] for proper interpolation!)
allowedz = [0.0001d0,0.0002d0,0.0005d0,0.0010d0,0.0020d0,$
            0.0040d0,0.0060d0,0.0080d0,0.0100d0,0.0140d0,$
            0.0170d0,0.0200d0,0.0300d0,0.0400d0,0.0600d0]
allowedy = [0.2490d0,0.2490d0,0.2490d0,0.2500d0,0.2520d0,$
            0.2560d0,0.2590d0,0.2630d0,0.2670d0,0.2730d0,$
            0.2790d0,0.2840d0,0.3020d0,0.3210d0,0.3560d0]
allowedx = 1d0-allowedy-allowedz
allowedmh = alog10((allowedz/allowedx)/zxsun)
allowedinitfeh = allowedmh ;; assumes alpha = 0

;; this must be in order for proper interpolation
sorted = sort(allowedinitfeh)
allowedinitfeh = allowedinitfeh[sorted]
allowedz = allowedz[sorted]
allowedy = allowedy[sorted]
allowedx = allowedx[sorted]

match = where(allowedinitfeh eq initfeh)
if match[0] eq -1 then begin
   print, '[Fe/H]_0 (' + strtrim(initfeh,2) + ') must be one of:'
   print, allowedinitfeh
   return,-1
endif

;; format the Z string
if round(allowedz[match[0]]*100d0) eq allowedz[match[0]]*100d0 then begin
   zstr = string(allowedz[match[0]],format='("Z",f0.2)')
endif else if round(allowedz[match[0]]*1000d0) eq allowedz[match[0]]*1000d0 then begin
   zstr = string(allowedz[match[0]],format='("Z",f0.3)')
endif else if round(allowedz[match[0]]*10000d0) eq allowedz[match[0]]*10000d0 then begin
   zstr = string(allowedz[match[0]],format='("Z",f0.4)')
endif
;; format the Y string
if round(allowedy[match[0]]*100d0) eq allowedy[match[0]]*100d0 then begin
   zstr += string(allowedy[match[0]],format='("Y",f0.2)')
endif else if round(allowedy[match[0]]*1000d0) eq allowedy[match[0]]*1000d0 then begin
   zstr += string(allowedy[match[0]],format='("Y",f0.3)')
endif

if 0 then begin
;; check for allowed alpha
allowedalpha = [0d0]
match = where(allowedalpha eq alpha)
if match[0] eq -1 then begin
   print, '[alpha/fe] (' + strtrim(alpha) + ') must be one of:'
   print, allowedalpha
   return,-1
endif
if alpha lt 0d0 then sign = 'm' $ 
else sign = 'p'
endif

;; Why does the mixing length change? This seems problematic...
if mstar le 0.7d0 then begin
   alphastr = 'A1.77_F7_' 
endif else begin
   alphastr = 'A1.74_F7_'
endelse

if 0 then begin
;; check for allowed vvcrit
allowedvvcrit = [0d0,0.4d0]
match = where(allowedvvcrit eq vvcrit)
if match[0] eq -1 then begin
   print, 'v/vcrit (' + strtrim(vvcrit,2) + ') must be one of:'
   print, allowedvvcrit
   return,-1
endif
vvcritstr = string(vvcrit,format='(f3.1)')
endif
vvcritstr = ''

;; mass grid points
allowedmass = [0.10d0,0.12d0,0.14d0,0.16d0,0.20d0,0.25d0,0.30d0,0.35d0,$
               0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,$
               0.80d0,0.85d0,0.90d0,0.95d0,1.00d0,1.05d0,1.10d0,$
               1.20d0,1.25d0,1.30d0,1.35d0,1.40d0,1.45d0,1.50d0,1.55d0,$
               1.60d0,1.65d0,1.70d0,1.75d0,1.80d0,1.85d0,1.95d0,2.00d0,$
               2.05d0,2.10d0,2.15d0,2.20d0,2.30d0,2.40d0,2.60d0,2.80d0,$
               3.00d0,3.20d0,3.40d0,3.60d0,3.80d0,4.00d0,4.20d0,4.40d0,$
;               4.60d0,4.80d0,5.00d0,5.20d0,5.40d0,       5.80d0,6.00d0,$
               4.60d0,4.80d0,5.00d0,5.20d0,5.40d0,5.60d0,5.80d0,6.00d0,$ ;; 5.6 is missing a grid point
;                      6.40d0,7.00d0,8.00d0,       10.0d0,       12.0d0,$ ;; 6.2, 9.0, 11 are missing a grid point
               6.20d0,6.40d0,7.00d0,8.00d0,9.00d0,10.0d0,11.0d0,12.0d0,$
;               14.0d0,16.0d0,18.0d0,20.0d0,       30.0d0,35.0d0,40.0d0,$ ;; 24 missing a grid point
               14.0d0,16.0d0,18.0d0,20.0d0,24.0d0,30.0d0,35.0d0,40.0d0,$

;                      50.0d0,55.0d0,60.0d0,       70.0d0,              $ ;; 45, 65, 75, 80 are missing grid points
              45.0d0,50.0d0,55.0d0,60.0d0,65.0d0,70.0d0,75.0d0,80.0d0,$
               90.0d0,100d0,120d0,150d0,200d0,250d0,300d0,350d0]

match = where(allowedmass eq mstar)
if match[0] eq -1 then begin
   print, 'mstar (' + strtrim(mstar,2) + ') must be one of:'
   print, allowedmass
   return,-1
endif
mstarstr = string(mstar,format='(f07.3)')

filename = zstr + 'OUT' + alphastr + 'M' + mstarstr + '.eep.idl'
eepfile = filepath(filename,$
                   root_dir=getenv('EXOFAST_PATH'),$
                   subdir=['parsec',zstr])

restore, eepfile
return, track

end
