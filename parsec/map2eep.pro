;; This reads in a PARSEC V1.2s track, resamples it onto the MIST
;; EEP grid, and stores it as an IDL save file for interpolation by
;; massradius_parsec.pro. It also generates diagnostic plots to
;; evaluate the quality of the EEP mapping

pro map2eep, trackname, display=display, overwrite=overwrite

;; most solar like track by default
if n_elements(trackname) eq 0 then trackname = 'Z0.014Y0.273/Z0.014Y0.273OUTA1.74_F7_M001.000.DAT'

;print, 'this is a problematic track for debugging! do not leave this here!'
;trackname = 'Z0.02Y0.284/Z0.02Y0.284OUTA1.77_F7_M000.090.DAT'


newtrackname = file_dirname(trackname) + path_sep() + file_basename(trackname,'.DAT') + '.eep.idl'
if ~keyword_set(overwrite) and file_test(newtrackname) then begin
   print, trackname + ' already done'
   return
endif

;; read in the tracks
readcol, trackname, model, mass, age, logl, logteff, logr, lograt, m_core_he,$
         m_core_c, h_cen, he_cen, c_cen, o_cen, lx, ly, lc, lneutr, lgrav, $
         h_sup, he_sup, c_sup, n_sup, o_sup, phase, skipline=1, $
         format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d'


constants = mkconstants()
eeps = dindgen(n_elements(model))
age /= 1d9 ;; Gyr
teff = 10^logteff ;; K
rstar = 10^logr/constants.rsun ;; R_sun

zxsun = 0.0207d0
y = he_sup
x = h_sup
z = 1d0 - x - y 
feh = alog10(z/x) - alog10(zxsun) ;; this is really [M/H] => alpha = 0

;; not enough sig figs! reject outliers...
;; x and y don't store enough sig figs
;; remove outliers for interpolation
;; smooth?
good = where(finite(feh))
djs_iterstat, feh[good], mean=mean, sigma=sigma, sigrej=5
bad = where(~finite(feh) or $
            feh gt (mean + 5d0*sigma) or $
            feh lt (mean - 5d0*sigma))
feh[bad] = !values.d_nan
good = where(finite(feh))
feh2 = interpol(feh[good], age[good], age) 
feh = feh2

;; determine the initial metalicity from Z and Y (encoded in the trackname)
initz = double((strsplit(file_basename(trackname),'ZY',/extract))[0])
inity = double((strsplit(file_basename(trackname),'ZY',/extract))[1])
initx = 1d0 - inity -initz
initfeh = alog10(initz/initx) - alog10(zxsun)

;; use MIST to find the corresponding EEPs for each point in the track
;; it doesn't have to be exact -- EEPs are fictional anyway
;; but the better this works, the more transferrable the EEP is
;; between the two models
tol = 1d-8
for i=0L, n_elements(model)-1 do begin
;   print, i+1, n_elements(model)

   mineep = 1d0
   maxeep = 4000d0

   while abs(maxeep-mineep) gt tol do begin
      eep = (mineep+maxeep)/2d0
      mistage = !values.d_nan

      ;; initfeh outside of MIST grid
      if initfeh gt 0.5d0 then initfeh = 0.5d0

      ;; mass outside of MIST grid
      if mass[i] lt 0.1d0 then mass[i] = 0.1d0
      if mass[i] gt 300d0 then mass[i] = 300d0

      chi2 = massradius_mist(eep, mass[i], initfeh, age[i], teff[i], rstar[i], feh[i], $
                             mistage=mistage, mistteff=mistteff, mistrstar=mistrstar,/allowold)

;      if finite(chi2) then begin
;         print, eep, mineep, maxeep, mistage, age[i], mistteff, teff[i], mistrstar, rstar[i], chi2
;         stop
;      endif

      if ~finite(mistage) then begin
         maxeep = eep
      endif else if age[i] lt mistage then begin
         maxeep = eep 
      endif else mineep = eep

   endwhile
   eeps[i] = eep

endfor

;spline=0
;quadratic=0
;lsquadratic=0
;; now interpolate the input track onto the MIST EEP grid
parsec_eeps = dindgen(ceil(max(eeps)))+1

;; trim all the repeated 1s, which screw up interpolation
bad = where(abs(eeps-1d0) lt tol,nbad)
if nbad gt 1 then begin
   eeps = eeps[bad[nbad-1]:*]
   age = age[bad[nbad-1]:*]
   logteff = logteff[bad[nbad-1]:*]
   feh = feh[bad[nbad-1]:*]
   logr = logr[bad[nbad-1]:*]
endif

if n_elements(parsec_eeps) ge 4 then begin
   eep_age = mcubic(age, eeps, parsec_eeps)
   eep_logteff = mcubic(logteff, eeps, parsec_eeps)
   eep_feh = mcubic(feh, eeps, parsec_eeps)
   eep_logrstar = mcubic(logr, eeps, parsec_eeps)
endif else if n_elements(parsec_eeps) eq 3 then begin
   eep_age = interpol(age, eeps, parsec_eeps,lsquadratic=lsquadradtic,quadratic=quadratic,spline=spline)
   eep_logteff = interpol(logteff, eeps, parsec_eeps,lsquadratic=lsquadradtic,quadratic=quadratic,spline=spline)
   eep_feh = interpol(feh, eeps, parsec_eeps,lsquadratic=lsquadradtic,quadratic=quadratic,spline=spline)
   eep_logrstar = interpol(logr, eeps, parsec_eeps,lsquadratic=lsquadradtic,quadratic=quadratic,spline=spline)
endif 

eep_teff = 10^eep_logteff
eep_rstar = 10^eep_logrstar/constants.rsun

set_plot, 'ps'
psname = file_dirname(trackname) + path_sep() + file_basename(trackname,'.DAT') + '.ps'
device, filename=psname, /color, bits=24
loadct, 39
red = 254
!p.font=0
!p.multi=[0,2,2]

;plot, eeps, age, xtitle='EEP', ytitle='Age',xrange=[0,max([eeps,parsec_eeps])], yrange=[min([age,eep_age]),max([age,eep_age])],/ys,/xs
;oplot, parsec_eeps, eep_age, color=red

plot, age, teff, xtitle='Age (Gyr)', ytitle='Teff',xrange=[0,max([age,eep_age])], yrange=[min([teff,eep_teff]),max([teff,eep_teff])],/ys,/xs
oplot, eep_age, eep_teff, color=red

plot, age, feh, xtitle='Age (Gyr)', ytitle='[Fe/H]',xrange=[0,max([age,eep_age])], yrange=[min([feh,eep_feh]),max([feh,eep_feh])],/ys,/xs
oplot, eep_age, eep_feh, color=red

plot, age, rstar, xtitle='Age (Gyr)', ytitle='Rstar',xrange=[0,max([age,eep_age])], yrange=[min([rstar,eep_rstar]),max([rstar,eep_rstar])],/ys,/xs,/ylog
oplot, eep_age, eep_rstar, color=red

logg = alog10(mass[0]/(rstar^2)*constants.gravitysun) ;; cgs
eep_logg = alog10(mass[0]/(eep_rstar^2)*constants.gravitysun) ;; cgs

plot, teff, logg, xtitle='Teff', ytitle='logg', xrange=[max(teff),min(teff)], yrange=[max(logg),min(logg)]
oplot, eep_teff, eep_logg, color=red

;plot, age, rstar, xtitle='Age',ytitle='Rstar'
;oplot, eep_age, eep_rstar, color=red

device, /close
set_plot, 'x'
!p.multi=0
if keyword_set(display) then spawn, 'gv ' + psname + ' &'

track = create_struct('age',eep_age,'teff',eep_teff,'ageweight',deriv(eep_age,parsec_eeps),'rstar',eep_rstar, 'feh', eep_feh)

save, track, filename=newtrackname

end
