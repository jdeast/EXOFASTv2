;; The SED constrains Teff, logg, [Fe/H], Extinction, and (Rstar/Distance)^2
function exofast_sed,fluxfile,teff,rstar,Av,d,logg=logg,met=met,alpha=alpha,f0=f0, fp0=fp0, ep0=ep0, verbose=verbose, psname=psname, pc=pc, rsun=rsun, interpfiles=interpfiles, logname=logname, fbol=fbol, debug=debug, redo=redo, flux=flux, oned=oned

common sed_block, klam, kkap, kapv, fp, ep, wp, widthhm, n, m, w1, kapp1

; fluxfile is the input set of fluxes to be fit
; a is the vector of input parameters
; f is the output model fluxes integrated over the bandpasses
; chi2 is output chi2

if not keyword_set(logg) then logg=4.5
if not keyword_set(met) then met=0.0
if not keyword_set(xmin) then xmin=0.1
if not keyword_set(xmax) then xmax=35.

;c=2.9979e14 ;; um/s
if n_elements(pc) eq 0 then pc=3.0857e18 ;; cm
if n_elements(rsun) eq 0 then rsun=6.96e10 ;; cm

;; don't do this every time
if n_elements(fp) eq 0 or keyword_set(redo) then begin
   mag2fluxconv,fluxfile,wp,widthhm,fp,ep,teff=teff 
   n=n_elements(wp)
   m=where(ep gt 0)
   
   readcol,filepath('extinction_law.ascii',root_dir=getenv('EXOFAST_PATH'),subdir='sed'),klam,kkap,/silent
   kapv = interpol(kkap,klam,0.55)

   w1 = findgen(24000)/1000+0.1 ;; wavelength scale on which to interpolate
   kapp1 = interpol(kkap,klam,w1)
endif

if keyword_set(oned) then begin
   ;; round to 0.5 dex in logg and met rather than interpolate
   logg1 = double(round(logg*2d0))/2d0
   exofast_interp_model,teff,logg1,met,w1,lamflam1temp,alpha=alpha,/next, interpfiles=interpfiles
endif else begin
   exofast_interp_model3d,teff,logg,met,w1,lamflam1temp,alpha=alpha,/next, interpfiles=interpfiles, logname=logname, verbose=verbose
endelse

if ~finite(lamflam1temp[0]) then begin
   f0 = !values.d_nan
   fp0 = fp
   ep0 = ep
   return, !values.d_infinity
endif

lamflam1=lamflam1temp*rstar*rstar*rsun*rsun/d/d/pc/pc
taul1 = kapp1/kapv/1.086*Av
extinct1 = exp(-taul1)

flux = lamflam1*extinct1

;fbol= total(flux)*(w1[1]-w1[0])
;print, 'sed derived: ' + strtrim(fbol,2)
;stop

f=fltarr(n)
for i=0,n-1 do begin
   ;; find the indices of the FWHM
   dummy = min((abs(w1-(wp[i]-widthhm[i]/2d0))),l)
   dummy = min((abs(w1-(wp[i]+widthhm[i]/2d0))),u)   
   f[i]=mean(flux[l:u])
endfor

if keyword_set(debug) or keyword_set(psname) eq 1 then begin
   ;; The Spectral Energy Distribution (black), with broad band
   ;; averages (blue circles) and broad band measurements (red). The error
   ;; bars in wavelength denote the bandwidth of the corresponding
   ;; filter and the error bars in flux denote the measurement
   ;; uncertainty.
   mydevice=!d.name
   if keyword_set(psname) then begin
      set_plot, 'PS'
      aspect_ratio=1.5
      xsize=10.5
      ysize=xsize/aspect_ratio
      !p.font=0
      device, filename=psname, /color, bits=24,/encapsulated
      device, xsize=xsize,ysize=ysize
      loadct, 39, /silent
      colors=[0,254,128,68]
      xtitle = exofast_textoidl('\lambda (\mum)')
      ytitle = exofast_textoidl('log \lambda F_\lambda (erg s^{-1} cm^{-2})')
      plotsym, 0, 0.5, /fill, color=colors[3]
   endif else begin
      device,window_state=win_state
      if win_state[5] eq 1 then wset, 5 $
      else window, 5, retain=2
      colors = ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x]
      xtitle = 'lambda (um)'
      ytitle = 'log(lambda F_lambda) (erg/s/cm^2)'
      plotsym, 0, 1.5, /fill, color=colors[3]
   endelse
   xmin = min(wp, max=xmax)
   xmax = 30
   xmin = 0.1
   ymin = alog10(min([f[m],fp[m]-ep[m]]))
   ymax = alog10(max([f[m],fp[m]+ep[m],smooth(flux,10)]))

   ;; overplot the model atmosphere
   plot, w1, alog10(smooth(flux,10)),/xlog,xtitle=xtitle,ytitle=ytitle, yrange=[ymin,ymax], xrange=[xmin,xmax], /xs;,/ys
   oplot, wp[m], alog10(f[m]), psym=8  

   ;; oploterror has too many dependencies
   for i=0, n_elements(m)-1 do begin
      ;; x error bar
      oplot, [wp[m[i]]-widthhm[m[i]]/2d0,wp[m[i]]+widthhm[m[i]]/2d0], alog10([fp[m[i]],fp[m[i]]]), color=colors[1]
      ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(wp[m[i]]-widthhm[m[i]]/2d0,alog10(fp[m[i]]),/to_device)
      xy2 = convert_coord(wp[m[i]]+widthhm[m[i]]/2d0,alog10(fp[m[i]]),/to_device)
      plots, [xy1[0],xy1[0]], [xy1[1]-ebw,xy1[1]+ebw], color=colors[1],/device
      plots, [xy2[0],xy2[0]], [xy2[1]-ebw,xy2[1]+ebw], color=colors[1],/device

      ;; y error bar
      oplot, [wp[m[i]],wp[m[i]]], alog10([fp[m[i]]-ep[m[i]],fp[m[i]]+ep[m[i]]]), color=colors[1]
      ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(wp[m[i]],alog10(fp[m[i]]-ep[m[i]]),/to_device)
      xy2 = convert_coord(wp[m[i]],alog10(fp[m[i]]+ep[m[i]]),/to_device)
      plots, [xy1[0]-ebw,xy1[0]+ebw], [xy1[1],xy1[1]], color=colors[1],/device
      plots, [xy2[0]-ebw,xy2[0]+ebw], [xy2[1],xy2[1]], color=colors[1],/device

   endfor

   if keyword_set(psname) then begin
      !p.font=-1
      !p.multi=0
      device, /close
      device, encapsulated=0

      residualfilename = file_dirname(psname) + path_sep() + file_basename(psname,'.eps') + '.residuals.txt'
      exofast_forprint, wp, f, fp, ep, f-ep, textout=residualfilename, comment='# Wavelength (um), model flux, flux, error, residuals (erg/s/cm^2)' 

   endif
   set_plot, mydevice


endif

f0 = f
fp0 = fp
ep0 = ep
;chi2 = total( (alog10(f[m]) - alog10(fp[m]))^2 / (ep[m]/fp[m]*alog(10)/2.5)^2)
chi2 = total(((f[m]-fp[m])/ep[m])^2)


return, chi2

end

