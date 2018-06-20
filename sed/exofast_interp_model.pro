;; t - teff
;; logg - the logg of the star
;; z  - metallicity
;; wav - wavelength, in ??
;; lamflam -
;; blamflam - 
;; dusty

;; marvels_interp_model, teff, logg, feh

pro exofast_interp_model,t,logg,feh,wav,lamflam,alpha=alpha,dusty=dusty,nextgen=nextgen,cond=cond,kurucz=kurucz,verbose=verbose, interpfiles=interpfiles

if t lt 800d0 or t gt 70000d0 then begin
   lamflam = !values.d_nan
   return
endif

;; loading files is **SLOW**
;; speed it up by keeping files loaded in this common block
common interp_model, model, specfiles, w1

if n_elements(alpha) eq 0 then alpha = 0.0

nextgen=1
if keyword_set(nextgen) then atdir=filepath('lte',root_dir=getenv("EXOFAST_PATH"),subdir=['sed','nextgenfin'])

if n_elements(model) eq 0 then begin
   specfiles = file_search(atdir + '*.spec.idl',count=nfiles)
   model = ptrarr(nfiles,/allocate_heap)
   w1 = findgen(24000)/1000+0.1 ;; wavelength scale on which to interpolate
endif

;; select the nearest models
allowedz = [-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.3,0.5]
junk = min(abs(allowedz-feh),ndx)
z = allowedz[ndx]

allowedalpha = [0d0,0.2d0,-0.2d0,0.4d0,0.6d0]
alphastr = ['+0.0','+0.2','-0.2','+0.4','+0.6']
nalpha = n_elements(allowedalpha)
junk = min(abs(allowedalpha-alpha),ndx)
a = allowedalpha[ndx]

allowedlogg = [-0.5d, 0d0, 0.5d0, 1d0, 1.5d0,  2d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0]
junk = min(abs(allowedlogg-logg),ndx)
g = allowedlogg[ndx]

allowedteff = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, $
               23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, $
               38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, $
               53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, $
               68, 69, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, $
               96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, $
               125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185,$
               190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,$
               310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430,$
               440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560,$
               570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690,$
               700]*100L
nteff = n_elements(allowedteff)
junk = min(abs(allowedteff-t),ndx)
if allowedteff[ndx] gt t then ndx -= 1
atnum1 = strtrim(allowedteff[ndx]/100,2)
atnum2 = strtrim(allowedteff[ndx+1]/100,2)

if n_params() lt 3 then begin
  print,'syntax: exofast_interp_model,t,logg,z,wav,lamflam,blamflam,dusty=dusty,nextgen=nextgen,cond=cond,kurucz=kurucz,verbose=verbose'
  retall
endif

if keyword_set(nextgen) then begin
   if g lt 0 then loggsign = '-'$
   else loggsign = '+'
   if z lt 0 then zsign = '-'$
   else zsign = '+'
   ;; ignore alpha -- most grid points only have one alpha anyway
endif

if keyword_set(nextgen) or keyword_set(dusty) or keyword_set(cond) then begin

   coolmatch = -1
   for i=0,-5,-1 do begin
      atnum = strtrim(allowedteff[ndx+i]/100,2)
      for j=0L, nalpha-1 do begin
         atend=string(loggsign, abs(g),format='(a,f3.1)')+string(zsign,abs(z),format='(a,f3.1)') + alphastr[j] +'.NextGen.spec.idl'
         coolname = atdir + atnum + atend
         coolmatch = (where(specfiles eq coolname))[0]
         if coolmatch ne -1 then break
      endfor
      if coolmatch ne -1 then break
   endfor

   hotmatch = -1
   for i=1,5 do begin
      atnum = strtrim(allowedteff[ndx+i]/100,2)
      for j=0L, nalpha-1 do begin
         atend=string(loggsign, abs(g),format='(a,f3.1)')+string(zsign,abs(z),format='(a,f3.1)') + alphastr[j] +'.NextGen.spec.idl'
         hotname = atdir + atnum + atend
         hotmatch = (where(specfiles eq hotname))[0]
         if hotmatch ne -1 then break
      endfor
      if hotmatch ne -1 then break
   endfor

   if hotmatch eq -1 or coolmatch eq -1 then begin
      lamflam = !values.d_nan
      if keyword_set(verbose) then message, "WARNING: No SED model for corresponding star (Teff=" + $
                                            string(t,logg,feh,alpha,format='(i5,", logg=",f5.2,", [Fe/H]=",f5.2,", alpha=",f5.2)') + '); rejecting model, which potentially imposes an unrealistic prior!',/continue
      return
   endif
endif
interpfiles = [coolname,hotname]

if keyword_set(nextgen) then begin
   ;; load the two models that bracket temperature of the requested
   ;; star for interpolation
   matches = [coolmatch,hotmatch]
   for i=0, 1 do begin

      ;; if it's never been read before, 
      ;; read it into a global (saved) variable
      if n_elements(*model[matches[i]]) eq 0 then begin

;         openr,myunit,files[i],/get_lun
;         readf,myunit,junk
;         readf,myunit,npoint
;         npoint=long(npoint)
;         lama1 = dblarr(npoint)
;         readf,myunit,lama1
;         lama1 = findgen(24000)/1000+0.1 ;; microns
;         lamflam1 = dblarr(npoint) ;; erg/s/cm^2
;         readf,myunit,lamflam1
;         close,myunit
;         free_lun,myunit

;; this interpolation has been moved to the model file generation step (see convertsed.pro)
;         lamflam1=flam1*lama1*1e-8 ;; convert from ?? to erg/s/cm^2
;         lama1=lama1*1e-4 ;; convert from angstroms to microns
;         lamflam1=interpol(lamflam1,lama1,w1)

         restore, specfiles[matches[i]]
         *model[matches[i]] = lamflam1 ;; populate the model variable for next time
      endif 
      if i eq 0 then lamflam1 = *model[matches[i]] $
      else lamflam2 = *model[matches[i]]
   endfor
endif

flams=dblarr(2,n_elements(lamflam1))
flams[0,*]=lamflam1
flams[1,*]=lamflam2

if keyword_set(nextgen) or keyword_set(dusty) or keyword_set(cond) then begin
   diff = fix(atnum2)-fix(atnum1)
   tfrac = 1d0-(fix((t+diff*100d0)/(diff*100d0))*diff*100d0 - t)/(diff*100d0)
endif

if keyword_set(kurucz) then begin
  tdiff=long(atnum2)-long(atnum1)
  tfrac=1.-(long((t+tdiff)/tdiff)*tdiff - t)/tdiff
endif

lamflam=reform(interpolate(flams,[tfrac],indgen(n_elements(lamflam1)),/grid))

wav=w1
   
if keyword_set(verbose) then begin 
   if keyword_set(ps) then begin
      loadct, 39, /silent
      colors=[0,254,68,128]
   endif else begin
      device,window_state=win_state
      if win_state[5] eq 1 then wset, 5 $
      else window, 5
      colors = ['ffffff'x,'ff0000'x,'00ff00'x,'0000ff'x]
   endelse
   
   ymax=max(lamflam1)
   plot,wav,smooth(lamflam1,15),/xlog,/ylog,yr=[ymax/20,ymax*1.2],xr=[0.2,30],/xs,/ys,xtitle=exofast_textoidl('\lambda (\mum)'), ytitle=exofast_textoidl('log \lambda F_\lambda (erg s^{-1} cm^{-2})')
   oplot,wav,smooth(lamflam2,15),col=colors[1]
   oplot,wav,smooth(lamflam,15),col=colors[2]

endif

end

