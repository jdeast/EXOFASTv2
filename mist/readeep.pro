function readeep, mstar, feh, vvcrit=vvcrit, alpha=alpha, rstar=rstar, teff=teff, age=age

if n_elements(alpha) eq 0 then alpha = 0d0
if n_elements(vvcrit) eq 0 then vvcrit = 0d0

;; check for allowed [Fe/H]
allowedfeh = [-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.25d0,$
              -1.0d0,-0.75d0,-0.5d0,-0.25d0,0d0,0.25d0,0.5d0]
match = where(allowedfeh eq feh)
if match[0] eq -1 then begin
   print, '[Fe/H] must be one of:'
   print, allowedfeh
   return,-1
endif
if feh lt 0d0 then sign = 'm' $ 
else sign = 'p'
fehstr = sign + string(abs(feh),format='(f4.2)')

;; check for allowed alpha
allowedalpha = [0d0]
match = where(allowedalpha eq alpha)
if match[0] eq -1 then begin
   print, '[alpha/fe] must be one of:'
   print, allowedalpha
   return,-1
endif
if alpha lt 0d0 then sign = 'm' $ 
else sign = 'p'
alphastr = sign + string(abs(alpha),format='(f3.1)')

;; check for allowed vvcrit
allowedvvcrit = [0d0,0.4d0]
match = where(allowedvvcrit eq vvcrit)
if match[0] eq -1 then begin
   print, 'v/vcrit must be one of:'
   print, allowedvvcrit
   return,-1
endif
vvcritstr = string(vvcrit,format='(f3.1)')

;; check for allowed mass
allowedmass = [0.10d0,0.15d0,0.20d0,0.25d0,0.30d0,0.35d0,0.40d0,$
               0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,0.80d0,0.85d0,$
               0.90d0,0.92d0,0.94d0,0.96d0,0.98d0,1.00d0,1.02d0,1.04d0,1.06d0,$
               1.08d0,1.10d0,1.12d0,1.14d0,1.16d0,1.18d0,1.20d0,1.22d0,1.24d0,$
               1.26d0,1.28d0,1.30d0,1.32d0,1.34d0,1.36d0,1.38d0,1.40d0,1.42d0,$
               1.44d0,1.46d0,1.48d0,1.50d0,1.52d0,1.54d0,1.56d0,1.58d0,1.60d0,$
               1.62d0,1.64d0,1.66d0,1.68d0,1.70d0,1.72d0,1.74d0,1.76d0,1.78d0,$
               1.80d0,1.82d0,1.84d0,1.86d0,1.88d0,1.90d0,1.92d0,1.94d0,1.96d0,$
               1.98d0,2.00d0,2.02d0,2.04d0,2.06d0,2.08d0,2.10d0,2.12d0,2.14d0,$
               2.16d0,2.18d0,2.20d0,2.22d0,2.24d0,2.26d0,2.28d0,2.30d0,2.32d0,$
               2.34d0,2.36d0,2.38d0,2.40d0,2.42d0,2.44d0,2.46d0,2.48d0,2.50d0,$
               2.52d0,2.54d0,2.56d0,2.58d0,2.60d0,2.62d0,2.64d0,2.66d0,2.68d0,$
               2.70d0,2.72d0,2.74d0,2.76d0,2.78d0,2.80d0,3.00d0,3.20d0,3.40d0,$
               3.60d0,3.80d0,4.00d0,4.20d0,4.40d0,4.60d0,4.80d0,5.00d0,5.20d0,$
               5.40d0,5.60d0,5.80d0,6.00d0,6.20d0,6.40d0,6.60d0,6.80d0,7.00d0,$
               7.20d0,7.40d0,7.60d0,7.80d0,8.00d0,9.00d0,10.00d0,11.00d0,$
               12.00d0,13.00d0,14.00d0,15.00d0,16.00d0,17.00d0,18.00d0,$
               19.00d0,20.00d0,22.00d0,24.00d0,26.00d0,28.00d0,30.00d0,$
               32.00d0,34.00d0,36.00d0,38.00d0,40.00d0,45.00d0,50.00d0,$
               55.00d0,60.00d0,65.00d0,70.00d0,75.00d0,80.00d0,85.00d0,$
               90.00d0,95.00d0,100.00d0,105.00d0,110.00d0,115.00d0,120.00d0,$
               125.00d0,130.00d0,135.00d0,140.00d0,145.00d0,150.00d0,175.00d0,$
               200.00d0,225.00d0,250.00d0,275.00d0,300.00d0]
match = where(allowedmass eq mstar)
if match[0] eq -1 then begin
   print, 'mstar must be one of:'
   print, allowedmass
   return,-1
endif
mstarstr = string(round(mstar*100d0),format='(i05)')

eepfile = filepath(mstarstr + 'M.track.eep.idl',$
                   root_dir=getenv('EXOFAST_PATH'),$
                   subdir=['mist','MIST_v1.0_tracks','feh_' + fehstr + '_afe_' + alphastr + $
                           '_vvcrit' + vvcritstr,'eeps'])

restore, eepfile
return, track

end
