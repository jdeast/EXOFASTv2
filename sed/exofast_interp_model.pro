;; t - teff
;; logg - the logg of the star
;; z  - metalicity
;; wav - wavelength, in ??
;; lamflam -
;; blamflam - 
;; dusty

;; marvels_interp_model, teff, logg, feh

pro exofast_interp_model,t,logg,feh,wav,lamflam,dusty=dusty,nextgen=nextgen,cond=cond,kurucz=kurucz,verbose=verbose

;; loading files is **SLOW**
;; speed it up by keeping files loaded in this common block
common interp_model, model, specfiles, w1

if n_elements(model) eq 0 then begin
   specfiles = file_search(getenv('EXOFAST_PATH') + '/sed/nextgenfin/*.spec',count=nfiles)
   model = ptrarr(nfiles,/allocate_heap)
   w1 = findgen(24000)/1000+0.1 ;; wavelength scale on which to interpolate
endif

if 0 then begin
yyz = [[0.00001d0,-3.29d0],$ 
     [0.00010d0,-2.29d0],$  
     [0.00040d0,-1.69d0],$  
     [0.00100d0,-1.29d0],$  
     [0.00400d0,-0.68d0],$ 
     [0.00700d0,-0.43d0],$  
     [0.01000d0,-0.27d0],$  
     [0.02000d0, 0.05d0],$  
     [0.04000d0, 0.39d0],$  
     [0.06000d0, 0.60d0],$  
     [0.08000d0, 0.78d0]]
z = interpol(yyz[0,*],yyz[1,*],feh)
endif
z = feh

nextgen=1

if n_params() lt 3 then begin
  print,'syntax: marvels_interp_model,t,logg,z,wav,lamflam,blamflam,dusty=dusty,nextgen=nextgen,cond=cond,kurucz=kurucz,verbose=verbose'
  retall
endif

atnum1 = strtrim(string(fix(t/100)), 2)
atnum2 = strtrim(string(fix((t+100)/100)), 2)

if keyword_set(nextgen) then atdir=getenv('EXOFAST_PATH') + '/sed/nextgenfin/lte'

if keyword_set(nextgen) then begin
  if z lt 0 then atend='-'+strn(logg,format='(f3.1)')+ strn(z,format='(f4.1)')+'.NextGen.spec'
  if z ge 0 then atend='-'+strn(logg,format='(f3.1)')+'-'+strn(z,format='(f3.1)')+'.NextGen.spec'
endif
if keyword_set(nextgen) then atend1=atend
if keyword_set(nextgen) then atend2=atend

if keyword_set(nextgen) or keyword_set(dusty) or keyword_set(cond) then begin
   while not file_test(atdir+atnum1+atend1) and atnum1 gt 0 do atnum1 = strtrim(string(atnum1-1),2)
   while not file_test(atdir+atnum2+atend2) and atnum2 lt 1000 do atnum2 = strtrim(string(atnum2+1),2)
   if not file_test(atdir+atnum2+atend2) or not file_test(atdir+atnum1+atend1) then begin
      lamflam = !values.d_nan
      return
   endif
endif

if keyword_set(nextgen) then begin

   ;; load the models on either side of the requested star
   match = (where(specfiles eq atdir+atnum1+atend1))[0]
   if n_elements(*model[match]) eq 0 then begin
      openr,myunit,atdir+atnum1+atend1,/get_lun
      readf,myunit,junk
      readf,myunit,npoint
      npoint=long(npoint)
      lama1 = dblarr(npoint)
      flam1 = dblarr(npoint)
      readf,myunit,lama1
      readf,myunit,flam1
      close,myunit
      free_lun,myunit
      lamflam1=flam1*lama1*1e-8
      lama1=lama1*1e-4
      lamflam1=interpol(lamflam1,lama1,w1)
      *model[match[0]] = lamflam1
   endif else lamflam1 = *model[match[0]]
   
   match = (where(specfiles eq atdir+atnum2+atend2))[0]
   if n_elements(*model[match]) eq 0 then begin
      openr,myunit,atdir+atnum2+atend2,/get_lun
      readf,myunit,junk
      readf,myunit,npoint
      npoint=long(npoint)
      lama2 = dblarr(npoint)
      flam2 = dblarr(npoint)
      readf,myunit,lama2
      readf,myunit,flam2
      close,myunit
      free_lun,myunit
      
      lamflam2=flam2*lama2*1e-8
      lama2=lama2*1e-4
      lamflam2=interpol(lamflam2,lama2,w1)
      *model[match[0]] = lamflam2
   endif else lamflam2 = *model[match[0]]   
endif

flams=dblarr(2,n_elements(lamflam1))
flams[0,*]=lamflam1
flams[1,*]=lamflam2

if keyword_set(nextgen) or keyword_set(dusty) or keyword_set(cond) then begin
   if fix(atnum2)-fix(atnum1) eq 1 then tfrac=1d0-(fix((t+100d0)/100d0)*100d0 - t)/100d0 $
   else if fix(atnum2)-fix(atnum1) eq 2 then tfrac=1d0-(fix((t+200d0)/200d0)*200d0 - t)/200d0 $
   else stop
endif

if keyword_set(kurucz) then begin
  tdiff=long(atnum2)-long(atnum1)
  tfrac=1.-(long((t+tdiff)/tdiff)*tdiff - t)/tdiff
endif

lamflam=reform(interpolate(flams,[tfrac],indgen(n_elements(lamflam1)),/grid))

wav=w1

if keyword_set(verbose) then begin 
   if keyword_set(ps) then begin
      set_plot, 'PS'
      device, filename='sed.ps',/color,bits=24
      loadct, 39, /silent
      colors=[0,254,68,128]
   endif else begin
      device,window_state=win_state
      if win_state[5] eq 1 then wset, 5 $
      else window, 5
      colors = ['ffffff'x,'ff0000'x,'00ff00'x,'0000ff'x]
   endelse
   
   ymax=max(lamflam1)
   plot,wav,smooth(lamflam1,15),/xlog,/ylog,yr=[ymax/20,ymax*1.2],xr=[0.2,30],/xs,/ys,xtitle=textoidl('\lambda (\mum)'), ytitle=textoidl('log \lambda F_\lambda (erg s^{-1} cm^{-2})')
   oplot,wav,smooth(lamflam2,15),col=colors[1]
   oplot,wav,smooth(lamflam,15),col=colors[2]
endif

end

