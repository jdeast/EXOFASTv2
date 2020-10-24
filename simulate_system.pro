pro simulate_system, minecc=minecc, maxecc=maxecc, ecc=ecc,$
                     minomega=minomega, maxomega=maxomega, omega=omega, $
                     mininitfeh=mininitfeh, maxinitfeh=maxinitfeh, initfeh=initfeh,$
                     mineep=mineep, maxeep=maxeep, eep=eep,$
                     minmstar=minmstar, maxmstar=maxmstar, mstar=mstar,$
                     minmp=minmp, maxmp=maxmp, mp=mp,$
                     minlogp=minlogp, maxlogp=maxlogp, logp=logp,$
                     mintc=mintc, maxtc=maxtc, tc=tc, $
                     mincosi=mincosi, maxcosi=maxcosi, cosi=cosi,$
                     nograzing=nograzing, neednottransit=neednottransit, i180=i180, $
                     mintime=mintime, maxtime=maxtime, cadence=cadence,$
                     maxiter=maxiter, prefix=prefix


;; define some physical constants
constants = mkconstants()
G = constants.GMSun/constants.RSun^3*constants.day^2 ;; R_sun^3/(m_sun*day^2)
Gmsun = constants.GMSun/constants.meter^3 ;; m^3/s^2
AU = constants.au/constants.rsun ;; R_sun
mjup = constants.gmjupiter/constants.gmsun ;; m_sun
rjup = constants.rjupiter/constants.rsun  ;; r_sun
mearth = constants.gmearth/constants.gmsun ;; m_sun
rearth = constants.rearth/constants.rsun  ;; r_sun
sigmaB = constants.sigmab/constants.lsun*constants.rsun^2 ;; Stefan-Boltzmann constant

if n_elements(dir) ne 1 then dir = './'
if n_elements(maxiter) ne 1 then maxiter = 50
niter = 0

restart:

niter ++
if niter gt maxiter then begin
   message, 'Drawn ' + strtrim(maxiter,2) + 'unphysical systems in a row. Your starting values/ranges may preclude a physical solution.'
endif

;; define ranges of parameters
if n_elements(minecc) ne 1 then minecc = 0d0
if n_elements(maxecc) ne 1 then maxecc = 1d0
if n_elements(ecc) ne 1 then ecc0 = minecc + (maxecc-minecc)*randomu(seed) $
else ecc0 = ecc

;; omega (radians)
if n_elements(minomega) ne 1 then minomega = 0d0
if n_elements(maxomega) ne 1 then maxomega = 2d0*!dpi
if n_elements(omega) ne 1 then omega0 = minomega + (maxomega-minomega)*randomu(seed) $
else omega0 = omega

;; mstar (solar masses)
if n_elements(minmstar) ne 1 then minmstar = 0.5d0
if n_elements(maxmstar) ne 1 then maxmstar = 2.0d0
if n_elements(mstar) ne 1 then mstar0 = minmstar + (maxmstar-minmstar)*randomu(seed) $
else mstar0 = mstar

;; Equal Evoluationary Phases
if n_elements(mineep) ne 1 then mineep = 202d0
if n_elements(maxeep) ne 1 then maxeep = 454d0
if n_elements(eep) ne 1 then eep0 = mineep + (maxeep-mineep)*randomu(seed) $
else eep0 = eep

;; metalicity at birth
if n_elements(mininitfeh) ne 1 then mininitfeh = -0.5d0
if n_elements(maxinitfeh) ne 1 then maxinitfeh = 0.5d0
if n_elements(initfeh) ne 1 then initfeh0 = mininitfeh + (maxinitfeh-mininitfeh)*randomu(seed) $
else initfeh0 = initfeh

;; mp, (mjup)
if n_elements(minmp) ne 1 then minmp = 0.001d0
if n_elements(maxmp) ne 1 then maxmp = 13d0
if n_elements(mp) ne 1 then mp0 = minmp +(maxmp-minmp)*randomu(seed) $
else mp0 = mp

;; period (days) -- log uniform
if n_elements(minperiod) ne 1 then minperiod = 0.3d0
if n_elements(maxperiod) ne 1 then maxperiod = 365d0
if n_elements(period) ne 1 then period0 = 10^(alog10(minperiod) + (alog10(maxperiod)-alog10(minperiod))*randomu(seed)) $
else period0 = period

;; span 1 year from TESS launch
if n_elements(mintime) ne 1 then mintime = julday(7,1,2018)
if n_elements(maxtime) ne 1 then maxtime = mintime + 365.25d0

;; tc -- draw uniformly within 1 period of the center of the allowed times
;; TODO: can go outside of user-supplied bounds! fix that
if n_elements(mintc) ne 1 then mintc = mintime
if n_elements(maxtc) ne 1 then maxtc = maxtime
if n_elements(tc) ne 1 then tc0 = (mintime+maxtime)/2d0 + (randomu(seed)-0.5d0)*period $
else t0 = tc

;; create a filename
if n_elements(filename) ne 1 then begin
   vcve = sqrt(1d0-ecc0^2)/(1d0+ecc0*sin(omega0))
   id = string(round(vcve*100d0),format='(i03)') + string(round(randomu(seed)*1d7),format='(i07)')
   path = dir + id + path_sep()
endif

;; derive Tp
phase=exofast_getphase(ecc0,omega0,/pri)  
tp = tc0 - period0*phase

;; derive Rstar, Teff, [Fe/H] based on MIST stellar evolutionary models
junk = massradius_mist(eep0, mstar0, initfeh0, 0.5d0, 6000d0, 1d0, 0d0, $
                       mistrstar=rstar, mistteff=teff, mistfeh=feh, mistage=age)

;; if the star is older than the universe, draw again
if age gt 13.86d0 then goto, restart

;; derive quadratic limb darkening parameters in TESS band based on
;; Claret tables and stellar values
if n_elements(band) ne 1 then band = 'TESS'
logg = alog10(mstar0/(rstar^2)*constants.gravitysun)
coeffs = quadld(logg, teff, feh, band)
u1 = coeffs[0]
u2 = coeffs[1]

;; derive Rplanet based on Chen & Kipping and Mplanet
mpsun = mp0*mjup                              ;; m_sun
rpsun = massradius_chen(mpsun/mearth)*rearth  ;; in r_sun
p = rpsun/rstar

;; derive semi-major axis based on Kepler's law
a = (period0^2*G*(mpsun+mstar0)/(4d0*!dpi^2))^(1d0/3d0) ;; r_sun
ar = a/rstar

;; not physical -- collides with star
if ecc0 gt 1d0 - (rstar+rpsun)/a then goto, restart

;; draw a random (transiting?) value for inclination
if n_elements(maxcosi) ne 1 then begin
   if keyword_set(neednottransit) then maxcosi = 1d0 $ ;; no restrictions on inclination
   else if keyword_set(nograzing) then maxcosi = (1d0-p)/ar*(1d0+ecc0*sin(omega0))/(1d0-ecc0^2) $ ;; must fully transit
   else maxcosi = (1d0+p)/ar*(1d0+ecc0*sin(omega0))/(1d0-ecc0^2) ;; must at least partially transit
endif
;; remove i = 90 +/- i degeneracy (not desired if fitting astrometry)
if n_elements(mincosi) ne 1 then begin
   if keyword_set(i180) then mincosi = -maxcosi $
   else mincosi = 0d0
endif
if n_elements(cosi) ne 1 then cosi0 = mincosi + (maxcosi-mincosi)*randomu(seed) $
else cosi0 = cosi

f0 = 1d0
variance = 0d0

inc = acos(cosi0)

;; transit duration
b = ar*cosi0*(1d0-ecc0^2)/(1d0+ecc0*sin(omega0))
t14 = period0/!dpi*asin(sqrt((1d0+p)^2 - b^2)/(sin(inc)*ar))*sqrt(1d0-ecc0^2)/(1d0+ecc0*sin(omega0))

;; derive RV semi-amplitude (m/s)
k = (2d0*!dpi*Gmsun/(period0*86400d0*(mstar0+mpsun)^2))^(1d0/3d0)*mpsun*sin(inc)/sqrt(1d0-ecc0^2) ;; m/s

file_mkdir, path

;; make a prior file
prefix = path + id + '.'
priorfile = prefix + 'priors'

urstar = rstar*0.035d0
uteff = teff*0.015d0
ufeh = 0.08d0

openw, lun, priorfile, /get_lun
printf, lun, rstar, urstar, format='("rstar ",f0.8,x,f0.8)' 
printf, lun, teff, uteff, format='("teff ",f0.8,x,f0.3)' 
printf, lun, feh, ufeh, format='("feh ",f0.8,x,f0.3)' 
printf, lun, initfeh0, format='("initfeh ",f0.8)' 
printf, lun, eep0, format='("eep ",f0.8)' 
printf, lun, mstar0, format='("mstar ",f0.8)' 
printf, lun, tc0, format='("tc ",f0.10)' 
printf, lun, omega0*180d0/!dpi, format='("omegadeg ",f0.8)' 
printf, lun, ecc0, format='("e ",f0.8)' 
printf, lun, mpsun, mpsun*0.05d0, format='("mpsun ",f0.8, x, f0.8)' 
printf, lun, p, format='("p ",f0.8)' 
printf, lun, cosi0, format='("cosi ",f0.8)' 
printf, lun, period0, format='("period ",f0.10)' 
printf, lun, u1, format='("u1 ",f0.8)' 
printf, lun, u2, format='("u2 ",f0.8)' 
printf, lun, variance, format='("variance ",f0.8)' 
printf, lun, f0, format='("f0 ",f0.8)'
free_lun, lun

;; create a tess 2-minute cadence LC
npoints = (maxtime-mintime)*24d0*30d0+1

;; this is the emitted time
emitted_time = mintime + dindgen(npoints)/(npoints-1)*365.25d0 ;; 1 year of 2 minute cadence continuously
q = mstar0/mpsun

;; convert times from target frame to observed (SSB) frame 
;; (they'll be converted back by exofast_tran)
observed_time = target2bjd(emitted_time, inc=inc, a=a,tp=tp, period=period0, e=ecc0, omega=omega0, q=q)

lc = exofast_tran(observed_time, inc, ar, tp, period0, ecc0, omega0, p, u1, u2, f0, rstar=rstar/au, q=q, au=au, c=c, tc=tc0, reflect=0d0, dilute=0d0, thermal=0d0)

;; add scatter to the LC                                                                                                           
if n_elements(err) eq 0 then err = 20d0/1d6
lc += randomn(seed,npoints)*err

;; save the full LC to a file 
filename = path + 'n20190101.TESS.TESS.dat'
forprint, observed_time, lc, replicate(err,npoints), format='(f0.6,x,f0.9,x,f0.9)', /nocomment, textout=filename
         
;; chop out the regions around the transit and save that
phase = (tesstime - tc) mod period
toohigh = where(phase gt period/2d0)
if toohigh[0] ne -1 then phase[toohigh] -= period
toolow = where(phase lt -period/2d0)
if toolow[0] ne -1 then phase[toolow] += period
keep = where(abs(phase) lt t14)
forprint, observed_time[keep], lc[keep], replicate(err,n_elements(keep)), format='(f0.6,x,f0.9,x,f0.9)', /nocomment, textout=filename+'.chopped'

;; TODO: make simulated RV/astrometric/SED/DT models...

end
