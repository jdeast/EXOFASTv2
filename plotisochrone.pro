pro plotisochrone, age, initfeh, teff=teff, rstar=rstar, mstar=mstar, feh=feh, parsec=parsec, plotmodel=plotmodel, epsname=epsname, isochrone_teff=isochrone_teff, isochrone_logg=isochrone_logg, isochrone_mstar=isochrone_mstar, isochrone_rstar=isochrone_rstar, xrange=xrange, yrange=yrange

if n_elements(xrange) ne 2 then xrange = [0,0]+!values.d_nan
if n_elements(yrange) ne 2 then yrange = [0,0]+!values.d_nan

;; default to the sun (this doesn't really matter)
if n_elements(teff) eq 0 then teff = 5788d0
if n_elements(rstar) eq 0 then rstar = 1d0
if n_elements(mstar) eq 0 then mstar = 1d0
if n_elements(feh) eq 0 then feh = 0d0
constants = mkconstants()
logg = alog10(mstar/(rstar^2)*constants.gravitysun)

nmstar = 1000d0
minmstar = alog10(0.1d0)
maxmstar = alog10(350d0)

isochrone_mstar = 10^(minmstar + dindgen(nmstar)/(nmstar-1)*(maxmstar-minmstar))
isochrone_rstar = dblarr(nmstar)+!values.d_nan
isochrone_teff = dblarr(nmstar)+!values.d_nan
age_isochrone = dblarr(nmstar)+!values.d_nan

tol = 1d-12
for i=0L, nmstar-1L do begin

   mineep = 1
   maxeep = 4000d0
   while abs(maxeep-mineep) gt tol do begin
      thiseep = (maxeep+mineep)/2d0
      
      model_rstar=!values.d_nan
      model_teff=!values.d_nan
      model_age=!values.d_nan
      
      if keyword_set(parsec) then begin
         junk = massradius_parsec(thiseep,isochrone_mstar[i],initfeh,age,teff[0],rstar[0],feh,$
                                  parsec_rstar=model_rstar,parsec_teff=model_teff,parsec_age=model_age,/allowold)
      endif else begin
         junk = massradius_mist(thiseep,isochrone_mstar[i],initfeh,age,teff[0],rstar[0],feh,$
                                mistrstar=model_rstar,mistteff=model_teff,mistage=model_age,/allowold)
      endelse
      
      if finite(junk) then begin
         if model_age gt age then maxeep = thiseep $
         else mineep = thiseep
         
;            if i ge 1 then $
;               print, massplot[i], mineep, maxeep, thiseep, parsec_age2, age, parsec_rstariso[i], parsec_teffiso[i], parsec_ageiso[i]
      endif else maxeep = thiseep
   endwhile
   
   if finite(junk) then begin
      if abs(age-model_age) lt 1d-6 then begin
         isochrone_rstar[i] = model_rstar
         isochrone_teff[i] = model_teff
         age_isochrone[i] = model_age
      endif
   endif         
endfor

isochrone_logg = alog10(isochrone_mstar/(isochrone_rstar^2)*constants.gravitysun)

ymin = max([isochrone_logg,logg],min=ymax) ;; plot range backwards
xmin = max([isochrone_teff,teff],min=xmax) ;; plot range backwards

;; plot ranges backwards
ymin = ceil(max(logg)*2d0)/2d0
ymax = floor(min(logg)*2d0)/2d0
xmin = ceil(max(teff)/1000d0)*1000d0
xmax = floor(min(teff)/1000d0)*1000d0

if finite(xrange[0]) then xmin=xrange[0]
if finite(xrange[1]) then xmax=xrange[1]
if finite(yrange[0]) then ymin=yrange[0]
if finite(yrange[1]) then ymax=yrange[1]

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
   window, 0, retain=2
   red = '0000ff'x
   symsize = 1
   xtitle='T_eff'
   ytitle='log g'
endelse



plot, isochrone_teff, isochrone_logg, xrange=[xmin,xmax], yrange=[ymin,ymax], xstyle=1, xticks=xticks, xminor=xminor, xtitle=xtitle, ytitle=ytitle

;; this overplots the input teff/logg and 
;; the nearest point on the isochrone
if keyword_set(plotmodel) then begin
   plotsym,0,/fill
   oplot, [teff], [logg], psym=8,symsize=0.5 ;; the model point
   for i=0L, n_elements(mstar)-1 do begin
      junk = min(abs(isochrone_mstar-mstar[i]),ndx)
      oplot, [isochrone_teff[ndx]], [isochrone_logg[ndx]], psym=2, color=red ;; overplot the time 
;      oplot, [isochrone_teff[ndx],teff[i]],[isochrone_logg[ndx],logg[i]],linestyle=1
   endfor
endif

if keyword_set(epsname) then begin
   !p.font=-1
   !p.multi=0
   device, /close
   device, encapsulated=0
endif

end
