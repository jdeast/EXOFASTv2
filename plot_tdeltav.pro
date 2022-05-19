pro plot_tdeltav, transitcsvfile

if n_elements(transitcsvfile) eq 0 then transitcsvfile = '/h/onion0/modeling/joey/new.fittransits.csv'
psname = file_basename(transitcsvfile,'.csv') + '.ps'
textname = file_basename(transitcsvfile,'.csv') + '.txt'

readcol, transitcsvfile, label, planet, epoch, tt, utt1, utt2, b, ub1, ub2, depth, udepth1, udepth2, p, up1, up2, delimiter=',',format='a,a,d,d,d,d,d,d,d,d,d,d,d,d,d'

readcol, filepath('filternames.txt', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist']), keivanname, mistname, format='a,a', comment='#',/silent

ntran = n_elements(label)
for i=0L, ntran-1 do begin

   prettyband = (strsplit(label[i],'()',/extract))[1]
   ;; translate pretty bands to machine-readable names
   if prettyband eq "u'" then band = 'SDSS_u' $
   else if prettyband eq "g'" then band = 'SDSS_g' $
   else if prettyband eq "r'" then band = 'SDSS_r' $
   else if prettyband eq "i'" then band = 'SDSS_i' $
   else if prettyband eq "z'" then band = 'SDSS_z' $
   else if prettyband eq "J" then band = '2MASS_J' $
   else if prettyband eq "H" then band = '2MASS_H' $
   else if prettyband eq "K" then band = '2MASS_Ks' $
   else if prettyband eq "Spit36" then band = 'IRAC_3.6' $
   else if prettyband eq "Spit45" then band = 'IRAC_4.5' $
   else if prettyband eq "Spit58" then band = 'IRAC_5.8' $
   else if prettyband eq "Spit80" then band = 'IRAC_8.0' $
   else band = prettyband

   filename = filepath(band + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
   if not file_test(filename) then begin
      match = (where(keivanname eq band))[0]
      filename = filepath(mistname[match] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
      if match[0] eq -1 or not file_test(filename) then begin
         print, band + ' not recognized. Bug in plot_tdelta?'
         stop
      endif        
   endif  
   restore, filename

   if i eq 0 then begin
      wp = filterproperties.lambda_eff/1d4
      widthhm = filterproperties.w_eff/1d4/2d0
      parr = p
      up1arr = up1
      up2arr = up2
   endif else begin
      wp = [wp,filterproperties.lambda_eff/1d4]
      widthhm = [widthhm,filterproperties.w_eff/1d4/2d0]
      parr = [parr,p]
      up1arr = [up1arr,up1]
      up2arr = [up2arr,up2]
   endelse

endfor


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


plotsym, 0,/fill
plot, wp, p, xrange=[min(wp-widthhm)*0.9,max(wp+widthhm)*1.1],yrange=[min(p-up2arr),max(p+up1arr)],xtitle=textoidl('Wavelength (\mum)'),ytitle=textoidl('R_P/R_*'),psym=8,symsize=0.5,/xlog,/xs

forprint, wp, p, widthhm, up2arr, up1arr, comment='#wavelength (um),rp/rs,bandpass half width (um), sigma_rp/rs_up, sigma_rp/rs_low',textout=textname

;; oploterror has too many dependencies; do it myself
for i=0, ntran-1 do begin
   ;; x error bar
   oplot, [wp[i]-widthhm[i],wp[i]+widthhm[i]], [p[i],p[i]], color=colors[1]
   ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i]-widthhm[i],p[i],/to_device)
   xy2 = convert_coord(wp[i]+widthhm[i],p[i],/to_device)
   plots, [xy1[0],xy1[0]], [xy1[1]-ebw,xy1[1]+ebw], color=colors[1],/device
   plots, [xy2[0],xy2[0]], [xy2[1]-ebw,xy2[1]+ebw], color=colors[1],/device
   
   ;; y error bar
   oplot, [wp[i],wp[i]], [parr[i]-up2arr[i],parr[i]+up1arr[i]], color=colors[1]
   ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i],parr[i]-up2arr[i],/to_device)
   xy2 = convert_coord(wp[i],parr[i]+up1arr[i],/to_device)
   plots, [xy1[0]-ebw,xy1[0]+ebw], [xy1[1],xy1[1]], color=colors[1],/device
   plots, [xy2[0]-ebw,xy2[0]+ebw], [xy2[1],xy2[1]], color=colors[1],/device
   
endfor

if keyword_set(psname) then begin
   device, /close
endif

end
