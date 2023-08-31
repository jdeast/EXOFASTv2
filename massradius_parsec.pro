;+
; NAME:
;   massradius_parsec
;
; PURPOSE: 
;   Interpolate the MIST stellar evolutionary models to derive Teff
;   and Rstar from mass, metallicity, and age. Intended to be a drop in
;   replacement for the Yonsie Yale model interpolation
;   (massradius_yy3.pro).
;
; CALLING SEQUENCE:
;   chi2 = massradius_parsec(mstar, feh, age, teff, rstar, $
;                          VVCRIT=vvcrit, ALPHA=alpha, SPAN=span,$
;                          PARSEC_RSTAR=mistrstar, PARSEC_TEFF=mistteff)
; INPUTS:
;
;    MSTAR  - The mass of the star, in m_sun
;    FEH    - The metallicity of the star [Fe/H]
;    AGE    - The age of the star, in Gyr
;    RSTAR  - The radius you expect; used to calculate a chi^2
;    TEFF   - The Teff you expect; used to calculate a chi^2
;    
; OPTIONAL INPUTS:
;   VVCRIT    - The rotational velocity normalized by the critical
;               rotation speed. Must be 0.0d0 or 0.4d0 (default 0.0d0).
;   ALPHA     - The alpha abundance. Must be 0.0 (default 0.0). A
;               placeholder for future improvements to MIST models.
;   SPAN      - The interpolation is done at the closest value +/-
;               SPAN grid points in the evolutionary tracks in mass,
;               age, metallicity. The larger this number, the longer it
;               takes. Default=1. Change with care.
;   EPSNAME   - A string specifying the name of postscript file to plot
;               the evolutionary track. If not specified, no plot is
;               generated.
;
; OPTIONAL KEYWORDS:
;   DEBUG     - If set, will plot the teff and rstar over the MIST
;               Isochrone.
;
; OPTIONAL OUTPUTS:
;   PARSEC_RSTAR - The rstar interpolated from the MIST models.
;   PARSEC_TEFF  - The Teff interpolated from the MIST models.
;
; RESULT:
;   The chi^2 penalty due to the departure from the MIST models,
;   assuming 3% errors in the MIST model values.
;
; COMMON BLOCKS:
;   MIST_BLOCK:
;     Loading EEPs (model tracks) is very slow. This common block
;     allows us to store the tracks in memory between calls. The first
;     call will take ~3 seconds. Subsequent calls that use the same
;     EEP files take 1 ms.
;
; EXAMPLE: 
;   ;; penalize a model for straying from the MIST models 
;   chi2 += massradius_parsec(mstar, feh, age, rstar=rstar, teff=teff)
;
; MODIFICATION HISTORY
; 
;  2018/01 -- Written, JDE
;-
function massradius_parsec, eep, mstar, initfeh, age, teff, rstar, feh, vvcrit=vvcrit, $
                            alpha=alpha, parsec_age=parsec_age, $
                            parsec_rstar=parsec_rstar, parsec_teff=parsec_teff,parsec_feh=parsec_feh,$
                            epsname=epsname, debug=debug,gravitysun=gravitysun,$
                            fitage=fitage, ageweight=ageweight, $
                            verbose=verbose, logname=logname, $
                            trackfile=trackfile, allowold=allowold, zxsun=zxsun,$
                            tefffloor=tefffloor, fehfloor=fehfloor, rstarfloor=rstarfloor, agefloor=agefloor

if n_elements(tefffloor) eq 0 then tefffloor = -1
if n_elements(fehfloor) eq 0 then fehfloor = -1
if n_elements(rstarfloor) eq 0 then rstarfloor = -1
if n_elements(agefloor) eq 0 then agefloor = -1

;massradius_parsec(398.26676d,1d0,0d0,4.593d0,5800d0,1d0,parsec_age=parsec_age, parsec_rstar=parsec_rstar, parsec_teff=parsec_teff) & print, parsec_age, parsec_teff, parsec_rstar

if n_elements(gravitysun) eq 0 then $                                           
   gravitysun = 27420.011d0 ;; cm/s^2                                           

;IDL> print, mstar, feh, age, teff, rstar, format='(f0.30)'
;mstar = 3.620807920892795017664411716396d0
;feh = -0.036250000000000004440892098501d0
;age = 4.758762521454782401519878476392d0
;teff = 6169.925000000000181898940354585648d0
;rstar = 1.416249999999999342747969421907d0

;; this common block allows us to store the tracks in memory between calls
common PARSEC_BLOCK, tracks, allowedmass, allowedinitfeh, nmass, nfeh, nvvcrit, nalpha

if mstar lt 0.1d0 or mstar gt 350d0 then begin
   if keyword_set(verbose) then printandlog, 'Mstar (' + strtrim(mstar,2) + ') is out of range [0.1,300]', logname
   return, !values.d_infinity
endif

if initfeh lt -4d0 or initfeh gt 0.69576805781432627d0 then begin
   if keyword_set(verbose) then printandlog, 'feh (' + strtrim(feh,2) + ') is out of range [-4,0.5]', logname
   return, !values.d_infinity
endif

if eep lt 1 then begin
   if keyword_set(verbose) then printandlog, 'EEP (' + strtrim(eep,2) + ') is out of range [1,infinity]', logname
   return, !values.d_infinity
endif

if n_elements(zxsun) eq 0 then zxsun=0.0207d0 ;; Bressan+ 2012

;; if this is the first call, initialize the tracks
if n_elements(tracks) eq 0 then begin

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

   ;; v/v_crit grid points
   allowedvvcrit = [0.0d0]
   
   ;; [alpha/Fe] grid points
   allowedalpha = [0d0]

   nmass = n_elements(allowedmass)
   nfeh = n_elements(allowedinitfeh)
   nvvcrit = n_elements(allowedvvcrit)
   nalpha = n_elements(allowedalpha)

   ;; each mass/metallicity points to a 3xN array for rstar and Teff for N ages
   ;; which we only load as needed (expensive)
   tracks = ptrarr(nmass, nfeh, nvvcrit, nalpha, /allocate_heap)

endif

;; load the closest track and one on either side for interpolation
junk = min(abs(mstar - allowedmass),massndx)
junk = min(abs(initfeh - allowedinitfeh),fehndx)
eepndx = round(eep)-1

;; don't interpolate over vvcrit or alpha
if n_elements(vvcrit) eq 0 then vvcritndx = 0 $
else begin
   vvcritndx = (where(allowedvvcrit eq vvcrit))[0]
   if vvcritndx eq -1 then message, 'vvcrit=' + strtrim(vvcrit,2) + ' not allowed'
endelse

if n_elements(alpha) eq 0 then alphandx = 0 $
else begin
   alphandx = (where(allowedalpha eq alpha))[0]
   if alphandx eq -1 then message, 'alpha=' + strtrim(alpha,2) + ' not allowed'
endelse

if n_elements(span) eq 0 then span = 1

;; make sure we don't overstep our bounds
if mstar lt allowedmass[massndx] then minmassndx = (massndx-1) > 0 $
else minmassndx = massndx < (nmass-2)

if initfeh lt allowedinitfeh[fehndx] then minfehndx = (fehndx-1) > 0 $
else minfehndx = fehndx < (nfeh-2)

if (eep-1) lt eepndx then mineepndx = (eepndx-1) > 0 $
else mineepndx = eepndx ;< (neep-2)

eepbox = mineepndx + [1,2]
mstarbox = allowedmass[minmassndx:minmassndx+1]
initfehbox = allowedinitfeh[minfehndx:minfehndx+1]

allfehs = dblarr(2,2,2)
allteffs = dblarr(2,2,2)
allages = dblarr(2,2,2)
allrstars = dblarr(2,2,2) 
allageweights = dblarr(2,2,2) 

for i=0, 1 do begin
   for j=0, 1 do begin
      
      ;; only load if it's not already loaded
      if n_elements(*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])) eq 0 then begin
         *(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]) = $ 
            readeep_parsec(allowedmass[i+minmassndx],allowedinitfeh[j+minfehndx], vvcrit=vvcrit, alpha=alpha, zxsun=zxsun)
      endif

      ;; extract the age, rstar, teff, and current feh of the matching track
      ages = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])).age
      rstars = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])).rstar
      teffs = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])).teff
      fehs = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])).feh
      ageweights = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx])).ageweight

      ;; make sure we don't overstep the EEP bounds
      neep = n_elements(ages)
      if eep gt neep then begin
         if keyword_set(verbose) then printandlog, 'EEP (' + strtrim(eep,2) + ') is out of range [1,' + strtrim(neep,2) + ']', logname
         return, !values.d_infinity
      endif

      mineepndx = floor(eep-1) < (neep-2)

      ;; populate the derived values
      allages[*,i,j] = ages[mineepndx:mineepndx+1]
      allrstars[*,i,j] = rstars[mineepndx:mineepndx+1]
      allteffs[*,i,j] = teffs[mineepndx:mineepndx+1]
      allfehs[*,i,j] = fehs[mineepndx:mineepndx+1]
      allageweights[*,i,j] = ageweights[mineepndx:mineepndx+1]
      
   endfor
endfor

x_eep = (eep-eepbox[0])/(eepbox[1]-eepbox[0])
y_mass = (mstar-mstarbox[0])/(mstarbox[1]-mstarbox[0])
z_feh = (initfeh-initfehbox[0])/(initfehbox[1]-initfehbox[0])

if float(!version.release) le 8.2 then begin
   parsec_age       = interpolate(allages  ,x_eep, y_mass, z_feh)
   parsec_rstar     = interpolate(allrstars,x_eep, y_mass, z_feh)
   parsec_teff      = interpolate(allteffs ,x_eep, y_mass, z_feh)
   parsec_feh       = interpolate(allfehs  ,x_eep, y_mass, z_feh)
   ageweight  = interpolate(allageweights  ,x_eep, y_mass, z_feh)
endif else begin
   parsec_age       = interpolate(allages  ,x_eep, y_mass, z_feh,/double)
   parsec_rstar     = interpolate(allrstars,x_eep, y_mass, z_feh,/double)
   parsec_teff      = interpolate(allteffs ,x_eep, y_mass, z_feh,/double)
   parsec_feh         = interpolate(allfehs  ,x_eep, y_mass, z_feh,/double)
   ageweight = interpolate(allageweights  ,x_eep, y_mass, z_feh,/double)
endelse

;; must be less than the age of the universe
if parsec_age lt 0 or (~keyword_set(allowold) and parsec_age gt 13.82d0) then begin
   if keyword_set(verbose) then printandlog, 'Age (' + strtrim(parsec_age,2) + ') is out of range',logname
   return, !values.d_infinity
endif

;; assume 3% model errors at 1 msun, 5% at 0.1 msun, 5% at 10 msun
percenterror = 0.03d0 + 0.02d0*alog10(mstar)^2

;; or overwrite with user supplied floors
if rstarfloor gt 0 then chi2 = ((parsec_rstar - rstar)/(rstarfloor*parsec_rstar))^2 $
else chi2 = ((parsec_rstar - rstar)/(percenterror*parsec_rstar))^2

if tefffloor gt 0 then chi2 += ((parsec_teff - teff)/(tefffloor*parsec_teff))^2 $
else chi2 += ((parsec_teff - teff)/(percenterror*parsec_teff))^2 

if fehfloor gt 0 then chi2 += ((parsec_feh - feh)/(fehfloor))^2 $
else chi2 += ((parsec_feh - feh)/(percenterror))^2 
              
if keyword_set(fitage) then begin
   if agefloor gt 0 then chi2 += ((parsec_age - age)/(agefloor*parsec_age))^2 $
   else chi2 += ((parsec_age - age)/(percenterror*parsec_age))^2
endif

;; plot it
if keyword_set(debug) or keyword_set(epsname) then begin
   mydevice=!d.name

   ;; prepare the plotting device
   if keyword_set(epsname) then begin
      ;; astrobetter.com tip on making pretty IDL plots
      set_plot, 'PS'
      aspect_ratio=1
      xsize=10.5/1.5
      ysize=xsize/aspect_ratio
      !p.font=0
      device, filename=epsname, /color, bits=24,/encapsulated
      device, xsize=xsize,ysize=ysize
      loadct, 39, /silent
      red = 254
      symsize = 0.33
      xtitle=exofast_textoidl('T_{eff}')
      ytitle=exofast_textoidl('log g_*')
   endif else begin
;      set_plot, 'X'
      red = '0000ff'x
      symsize = 1
      ;device,window_state=win_state
      ;if win_state[29] eq 1 then wset, 29 $
      ;else window, 29, retain=2
      xtitle='T_eff'
      ytitle='log g'
   endelse

   ;; interpolate the entire track to plot the mass track
   mineep = min([202,eep])
   maxeep = 808
   neep = maxeep - mineep + 1
   eepplot = mineep + dindgen(neep)
   parsec_rstariso = dblarr(neep)
   parsec_teffiso = dblarr(neep)
   parsec_ageiso = dblarr(neep)

   for i=0, neep-1 do begin
      junk = massradius_parsec(eepplot[i],mstar,initfeh,age,teff,rstar,feh,parsec_rstar=parsec_rstar2,parsec_teff=parsec_teff2,parsec_age=parsec_age2)
      parsec_rstariso[i] = parsec_rstar2
      parsec_teffiso[i] = parsec_teff2
      parsec_ageiso[i] = parsec_age2
   endfor

;; this is to plot an isochrone rather than a mass track
if 0 then begin
   maxlogmass = alog10(10d0)
   minlogmass = alog10(0.1d0)
   nmass = 10000d0
   massplot = 10^(minlogmass + dindgen(nmass)/(nmass-1d0)*(maxlogmass-minlogmass))
   tol = 1d-8
   parsec_rstariso = dblarr(nmass)+!values.d_nan
   parsec_teffiso = dblarr(nmass)+!values.d_nan
   parsec_ageiso = dblarr(nmass)+!values.d_nan
   for i=0, nmass-1 do begin
      mineep = 0
      maxeep = 1808
      while abs(maxeep-mineep) gt tol do begin
         thiseep = (maxeep+mineep)/2d0
         
         parsec_rstar=!values.d_nan
         parsec_teff=!values.d_nan
         parsec_age=!values.d_nan

         junk = massradius_parsec(thiseep,massplot[i],initfeh,age,teff,rstar,feh,parsec_rstar=parsec_rstar2,parsec_teff=parsec_teff2,parsec_age=parsec_age2,/allowold)

         if finite(junk) then begin
            if parsec_age2 gt age then maxeep = thiseep $
            else mineep = thiseep
            
;            if i ge 1 then $
;               print, massplot[i], mineep, maxeep, thiseep, parsec_age2, age, parsec_rstariso[i], parsec_teffiso[i], parsec_ageiso[i]
         endif else maxeep = eep
;         wait, 0.05
      endwhile

      if finite(junk) then begin
         if abs(age-parsec_age2) lt 1d-6 then begin
            parsec_rstariso[i] = parsec_rstar2
            parsec_teffiso[i] = parsec_teff2
            parsec_ageiso[i] = parsec_age2
;      parsec_eep[i] = thiseep
         endif
      endif

;      print, massplot[i], mineep, maxeep, thiseep, parsec_age2, age, parsec_rstariso[i], parsec_teffiso[i], parsec_ageiso[i]
;stop
         
   endfor
endif


   if n_elements(trackfile) ne 0 then $
      exofast_forprint, parsec_teffiso, parsec_rstariso, parsec_ageiso, format='(f0.5,x,f0.5,x,f0.5)', comment='#teff, rstar, age', textout=trackfile

   ;; make a publication-ready plot of the PARSEC track -- 
   ;; Teff vs logg at a fixed mass over time
   loggplottrack =  alog10(mstar/(parsec_rstariso^2)*gravitysun)
   teffplottrack = parsec_teffiso
   loggplot =  alog10(mstar/(rstar^2)*gravitysun)

   use = where(loggplottrack gt 3 and loggplottrack lt 5); and eepplot ge min([202,eep]))
   if use[0] eq -1 then use = lindgen(n_elements(loggplottrack))

   xmin=max([teffplottrack[use],teff],min=xmax) ;; plot range backwards

   ;; increase xrange so there are 4 equally spaced ticks that land on 100s
   xticks = 3
   xmin = ceil(xmin/100)*100
   xmax = floor(xmax/100)*100
   repeat begin
      spacing = ceil((xmin-xmax)/3d0/100d0)*100d0
      if (xmin-xmax)/spacing ne xticks then begin
         xmin += 100d0
         xmax -= 100d0
      endif
   endrep until (xmin-xmax)/spacing eq xticks
   xminor = spacing/100d0
  
   ymax = min([loggplot,3,5,loggplottrack[use]],max=ymin) ;; plot range backwards
   plot, teffplottrack[use], loggplottrack[use],xtitle=xtitle,ytitle=ytitle, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, xticks=xticks, xminor=xminor
   plotsym,0,/fill
   oplot, [teff], [loggplot], psym=8,symsize=0.5 ;; the model point
   junk = min(abs(eepplot-eep),ndx)
   oplot, [teffplottrack[ndx]],[loggplottrack[ndx]], psym=2, color=red ;; overplot the time 
   if keyword_set(epsname) then begin
      !p.font=-1
      !p.multi=0
      device, /close
      device, encapsulated=0
   endif

   set_plot, mydevice

endif

return, chi2

end
