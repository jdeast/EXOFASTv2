pro plotsed, 

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
   ;; black, red, green, blue, purple, orange, yellow 
   colors=[0,254,128,68,32,208,192] 
   xtitle = exofast_textoidl('\lambda (\mum)')
   ytitle = exofast_textoidl('log \lambda F_\lambda (erg s^{-1} cm^{-2})')
   ytitle2 = exofast_textoidl('Res \sigma')
   
   symsize = 0.5
   plotsym, 0, symsize, /fill, color=colors[3]
   
   position1 = [0.23, 0.40, 0.95, 0.95]    ;; data plot
   position2 = [0.23, 0.20, 0.95, 0.40]    ;; residual plot
   
endif else begin
   device,window_state=win_state
   if win_state[5] eq 1 then wset, 5 $
   else window, 5, retain=2
   colors = ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x]
   xtitle = 'lambda (um)'
   ytitle = 'log(lambda F_lambda) (erg/s/cm^2)'
   symsize=1.5
   plotsym, 0, symsize, /fill, color=colors[3]
   
   position1 = [0.07, 0.22, 0.97, 0.95]    ;; data plot
   position2 = [0.07, 0.07, 0.97, 0.22]    ;; residual plot
endelse
ncolors = n_elements(colors)

wp = transpose(filterprops.lambda_eff)/1d4        ;; effective wavelength, in um
widthhm = transpose(filterprops.w_eff)/1d4/2d0    ;; half width, in um

;; convert everything to AB mag scale
vega = where(filterprops.type eq 'Vega')
abmags = mags
modelabmags = blendmag          ;modelmags
if vega[0] ne -1 then begin
   abmags[vega] += filterprops[vega].vegaab
   modelabmags[vega] += filterprops[vega].vegaab
endif

;; compute the absolute fluxes
zp = 3631d0*3d-9/wp
flux = zp*10^(-0.4d0*abmags)
fluxerr = flux*alog(10d0)/2.5d0*errs
modelflux = zp*10^(-0.4d0*modelabmags)

;; define the limits; plot the model flux
xmin = min(wp, max=xmax)
xmax = 30
xmin = 0.1
ymin = alog10(min([modelflux,flux-fluxerr])) ;,reform(atmospheres,n_elements(atmospheres))]))
ymax = alog10(max([modelflux,flux+fluxerr,total(atmospheres,1)]))
if finite(range[0]) then xmin = range[0]
if finite(range[1]) then xmax = range[1]
if finite(range[2]) then ymin = range[2]
if finite(range[3]) then ymax = range[3]

plot, [0], [0], /xlog, ytitle=ytitle, yrange=[ymin,ymax], xrange=[xmin,xmax], /xs, position=position1, xtickformat='(A1)'

;; plot each individual star
for i=0L, nstars-1 do begin
   oplot, wavelength, alog10(smooth(atmospheres[i,*],10)), color=colors[(i+1) mod ncolors]
endfor

;; set colors for each star and all stars
pointcolors = dblarr(nbands)
for i=0L, nbands-1 do begin
   if total(blend[i,*]) eq 1d0 then pointcolors[i] = colors[((where(blend[i,*] eq 1L))[0]+1) mod ncolors]
   if total(blend[i,*]) eq nstars then pointcolors[i] = colors[0]
endfor

;; plot each combination of blended stars 
;; except 1 and all, handled separately
starnames = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
legendlabels = starnames[lindgen(nstars)]
legendcolors = colors[(lindgen(nstars)+1) mod ncolors]
colorndx = nstars+1
for i=0L, nbands-1L do begin
   if total(blend[i,*]) gt 1 and total(blend[i,*]) ne nstars then begin
      for j=0L, i-1 do begin
         ;; if no other star before has this combination
         match = where(legendlabels eq strjoin(starnames[where(blend[i,*])],'+'),nmatch)
         if nmatch eq 0 then begin
            ;; plot it
            blended_atmosphere = transpose(atmospheres[0,*])*0d0
            for k=0L, nstars-1 do begin
               if blend[i,k] then blended_atmosphere += atmospheres[k,*]
            endfor
            color = colors[colorndx mod ncolors]
            oplot, wavelength, alog10(smooth(blended_atmosphere,10)), color=colors[colorndx mod ncolors]
            legendlabels = [legendlabels,strjoin(starnames[where(blend[i,*])],'+')]
            legendcolors = [legendcolors,color]
            colorndx++
            pointcolors[i] = color
         endif
      endfor
   endif
endfor
;; plot all blended together
oplot, wavelength, alog10(smooth(total(atmospheres,1),10)),color=colors[0]
legendlabels = [legendlabels,strjoin(starnames[0:nstars-1],'+')]
legendcolors = [legendcolors,colors[0]]
if nstars gt 1 then exofast_legend, legendlabels, textcolors=legendcolors,/top,/right,charsize=0.75

;; plot bands
oplot, wp, alog10(modelflux), psym=8

;; oploterror has too many dependencies; do it myself
for i=0, nbands-1 do begin
   ;; x error bar
   oplot, [wp[i]-widthhm[i],wp[i]+widthhm[i]], alog10([flux[i],flux[i]]), color=colors[1]
   ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i]-widthhm[i],alog10(flux[i]),/to_device)
   xy2 = convert_coord(wp[i]+widthhm[i],alog10(flux[i]),/to_device)
   plots, [xy1[0],xy1[0]], [xy1[1]-ebw,xy1[1]+ebw], color=colors[1],/device
   plots, [xy2[0],xy2[0]], [xy2[1]-ebw,xy2[1]+ebw], color=colors[1],/device
   
   ;; y error bar
   oplot, [wp[i],wp[i]], alog10([flux[i]-fluxerr[i],flux[i]+fluxerr[i]]), color=colors[1]
   ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i],alog10(flux[i]-fluxerr[i]),/to_device)
   xy2 = convert_coord(wp[i],alog10(flux[i]+fluxerr[i]),/to_device)
   plots, [xy1[0]-ebw,xy1[0]+ebw], [xy1[1],xy1[1]], color=colors[1],/device
   plots, [xy2[0]-ebw,xy2[0]+ebw], [xy2[1],xy2[1]], color=colors[1],/device
   
endfor

;; now plot the residuals below
residuals = (flux-modelflux)/fluxerr
res_errhi = (flux-modelflux+fluxerr)/fluxerr
res_errlo = (flux-modelflux-fluxerr)/fluxerr

;; round to 0.5
ymin = floor(min(res_errlo)/0.5)*0.5
ymax = ceil(max(res_errhi)/0.5)*0.5

;; make yrange symmetric
if abs(ymin) gt abs(ymax) then ymax =  abs(ymin)
if abs(ymax) gt abs(ymin) then ymin = -abs(ymax)

if finite(range[4]) then ymin = range[4]
if finite(range[5]) then ymax = range[5]

plot, [0],[0], position=position2, /noerase, $
      xrange=[xmin,xmax], xtitle=xtitle, /xlog, $
      yrange=[ymin,ymax]/0.7, ytitle=ytitle2,$ ;ytitle='O-C',$;ytitle=textoidl('O-C (\sigma)'), $
      /xstyle, /ystyle, yminor=2,yticks=2,ytickv=[ymin,0d0,ymax]
oplot, [xmin,xmax], [0d0,0d0], linestyle=1,color=colors[1]
for i=0L, nbands-1 do begin
   plotsym, 0, symsize, /fill, color=pointcolors[i]
   oplot, [wp[i]], [residuals[i]], psym=8 ;, color=pointcolors[i]
   
   ;; x error bar
   oplot, [wp[i]-widthhm[i],wp[i]+widthhm[i]], [residuals[i],residuals[i]], color=pointcolors[i]
   ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i]-widthhm[i],residuals[i],/to_device)
   xy2 = convert_coord(wp[i]+widthhm[i],residuals[i],/to_device)
   
   if xy1[0] lt 0 then xy1[0] = 0d0
   if xy1[0] gt !d.x_size then xy1[0] = !d.x_size
   if xy1[1] lt 0 then xy1[1] = 0d0
   if xy1[1] gt !d.y_size then xy1[1] = !d.y_size
   if xy2[0] lt 0 then xy2[0] = 0d0
   if xy2[0] gt !d.x_size then xy2[0] = !d.x_size
   if xy2[1] lt 0 then xy2[1] = 0d0
   if xy2[1] gt !d.y_size then xy2[1] = !d.y_size
   
   plots, [xy1[0],xy1[0]], [xy1[1]-ebw,xy1[1]+ebw], color=pointcolors[i],/device
   plots, [xy2[0],xy2[0]], [xy2[1]-ebw,xy2[1]+ebw], color=pointcolors[i],/device
   
   ;; y error bar
   oplot, [wp[i],wp[i]], [res_errlo[i],res_errhi[i]], color=pointcolors[i]
   ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
   xy1 = convert_coord(wp[i],res_errlo[i],/to_device)
   xy2 = convert_coord(wp[i],res_errhi[i],/to_device)
   plots, [xy1[0]-ebw,xy1[0]+ebw], [xy1[1],xy1[1]], color=pointcolors[i],/device
   plots, [xy2[0]-ebw,xy2[0]+ebw], [xy2[1],xy2[1]], color=pointcolors[i],/device
   
endfor


;; clean up the postscript device
if keyword_set(psname) then begin
   !p.font=-1
   !p.multi=0
   device, /close
   device, encapsulated=0
   
   ;; create a residual file
   residualfilename = file_dirname(psname) + path_sep() + 'modelfiles' + path_sep() + file_basename(psname,'.eps') + '.residuals.txt'
   
   startxt = strarr(nbands)
   for i=0L, nbands-1 do startxt[i] = strjoin(strtrim(where(blend[i,*]),2),',')
   
   exofast_forprint, filterprops.name, wp, widthhm, flux, fluxerr, modelflux, flux-modelflux,startxt, textout=residualfilename, $
                     comment='# Filtername, Center wavelength (um), half bandpass (um), flux, error, modelflux, residuals (erg/s/cm^2), star indices', $
                     format='(a20,x,f0.6,x,f0.6,x,e0.6,x,e0.6,x,e0.6,x,e0.6,x,a)'
   
endif
set_plot, mydevice

end
