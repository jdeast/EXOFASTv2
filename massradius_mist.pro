;+
; NAME:
;   massradius_mist
;
; PURPOSE: 
;   Interpolate the MIST stellar evolutionary models to derive Teff
;   and Rstar from mass, metalicity, and age. Intended to be a drop in
;   replacement for the Yonsie Yale model interpolation
;   (massradius_yy3.pro).
;
; CALLING SEQUENCE:
;   chi2 = massradius_mist(mstar, feh, age, teff, rstar, $
;                          VVCRIT=vvcrit, ALPHA=alpha, SPAN=span,$
;                          MISTRSTAR=mistrstar, MISTTEFF=mistteff)
; INPUTS:
;
;    MSTAR  - The mass of the star, in m_sun
;    FEH    - The metalicity of the star [Fe/H]
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
;               age, metalicity. The larger this number, the longer it
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
;   MISTRSTAR - The rstar interpolated from the MIST models.
;   MISTTEFF  - The Teff interpolated from the MIST models.
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
;   chi2 += massradius_mist(mstar, feh, age, rstar=rstar, teff=teff)
;
; MODIFICATION HISTORY
; 
;  2018/01 -- Written, JDE
;-
function massradius_mist, eep, mstar, initfeh, age, teff, rstar, feh, vvcrit=vvcrit, alpha=alpha, $
                               mistage=mistage, mistrstar=mistrstar, mistteff=mistteff, mistfeh=mistfeh,$
                               epsname=epsname, debug=debug, gravitysun=gravitysun, fitage=fitage, ageweight=ageweight

if n_elements(gravitysun) eq 0 then $                                           
   gravitysun = 27420.011d0 ;; cm/s^2                                           

;IDL> print, mstar, feh, age, teff, rstar, format='(f0.30)'
;mstar = 3.620807920892795017664411716396d0
;feh = -0.036250000000000004440892098501d0
;age = 4.758762521454782401519878476392d0
;teff = 6169.925000000000181898940354585648d0
;rstar = 1.416249999999999342747969421907d0


;; this common block allows us to store the tracks in memory between calls
common mist_block, tracks, allowedmass, allowedinitfeh, nmass, nfeh, nvvcrit, nalpha

if mstar lt 0.1d0 or mstar gt 300d0 then return, !values.d_infinity
if initfeh lt -4d0 or initfeh gt 0.5d0 then return, !values.d_infinity
if eep lt 1 then return, !values.d_infinity

;; if this is the first call, initialize the tracks
if n_elements(tracks) eq 0 then begin

   ;; mass grid points
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

   ;; [Fe/H] grid points
   allowedinitfeh = [-4.0d0,-3.5d0,-3.0d0,-2.5d0,-2.0d0,-1.5d0,-1.25d0,$
                     -1.0d0,-0.75d0,-0.5d0,-0.25d0,0d0,0.25d0,0.5d0]

   ;; v/v_crit grid points
   allowedvvcrit = [0.0d0, 0.4d0]
   
   ;; [alpha/Fe] grid points
   allowedalpha = [0d0]

   nmass = n_elements(allowedmass)
   nfeh = n_elements(allowedinitfeh)
   nvvcrit = n_elements(allowedvvcrit)
   nalpha = n_elements(allowedalpha)

   ;; each mass/metalicity points to a 3xN array for rstar and Teff for N ages
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
            readeep(allowedmass[i+minmassndx],allowedinitfeh[j+minfehndx], vvcrit=vvcrit, alpha=alpha)
      endif

      ;; extract the age, rstar, teff, and current feh of the matching track
      ages = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]))[0,*]
      rstars = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]))[1,*]
      teffs = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]))[2,*]
      fehs = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]))[3,*]
      ageweights = (*(tracks[i+minmassndx,j+minfehndx,vvcritndx,alphandx]))[4,*]

      ;; make sure we don't overstep the EEP bounds
      neep = n_elements(ages)
      if eep gt neep then return, !values.d_infinity
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
   mistage       = interpolate(allages  ,x_eep, y_mass, z_feh)
   mistrstar     = interpolate(allrstars,x_eep, y_mass, z_feh)
   mistteff      = interpolate(allteffs ,x_eep, y_mass, z_feh)
   mistfeh       = interpolate(allfehs  ,x_eep, y_mass, z_feh)
   ageweight = interpolate(allageweights  ,x_eep, y_mass, z_feh)
endif else begin
   mistage       = interpolate(allages  ,x_eep, y_mass, z_feh,/double)
   mistrstar     = interpolate(allrstars,x_eep, y_mass, z_feh,/double)
   mistteff      = interpolate(allteffs ,x_eep, y_mass, z_feh,/double)
   mistfeh       = interpolate(allfehs  ,x_eep, y_mass, z_feh,/double)
   ageweight = interpolate(allageweights  ,x_eep, y_mass, z_feh,/double)
endelse

;; must be less than the age of the universe
if mistage lt 0 or mistage gt 13.82d0 then return, !values.d_infinity

;; assume 3% model errors
chi2 = ((mistrstar - rstar)/(0.03d0*mistrstar))^2
chi2 += ((mistteff - teff)/(0.03d0*mistteff))^2
if keyword_set(fitage) then chi2 += ((mistage - age)/(0.03d0*mistage))^2
chi2 += ((mistfeh - feh)/(0.03d0))^2

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
      red = '0000ff'x
      symsize = 1
      device,window_state=win_state
      if win_state[29] eq 1 then wset, 29 $
      else window, 29, retain=2
      xtitle='T_eff'
      ytitle='log g'
   endelse

   ;; interpolate the entire track to plot the isochrone
   mineep = 150
   maxeep = 808
   neep = maxeep - mineep + 1
   eepplot = mineep + dindgen(neep)
   mistrstariso = dblarr(neep)
   mistteffiso = dblarr(neep)
   for i=0, neep-1 do begin
      junk = massradius_mist(eepplot[i],mstar,initfeh,age,teff,rstar,feh,mistrstar=mistrstar,mistteff=mistteff)
      mistrstariso[i] = mistrstar
      mistteffiso[i] = mistteff
   endfor

   ;; make a publication-ready plot of the YY track -- Teff vs logg
   loggplottrack =  alog10(mstar/(mistrstariso^2)*gravitysun)
   teffplottrack = mistteffiso
   loggplot =  alog10(mstar/(rstar^2)*gravitysun)

   use = where(loggplottrack gt 3 and loggplottrack lt 5)
   xmin=max(teffplottrack[use],min=xmax) ;; plot range backwards
   ymin = min([loggplot,3,5],max=ymax)
   
   plot, teffplottrack, loggplottrack,xtitle=xtitle,ytitle=ytitle, xrange=[xmin,xmax], yrange=[ymin,ymax];, xtickinterval=1500
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
