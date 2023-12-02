;; The SED constrains Teff, logg, [Fe/H], Extinction, and (Rstar/Distance)^2
function exofast_multised,teff, logg, feh, av, distance, lstar, errscale, sedfile, alpha=alpha, debug=debug, psname=psname, range=range, specphotpath=specphotpath, logname=logname,redo=redo,blend0=blend0,rstar=rstar, sperrscale=sperrscale


nstars = n_elements(teff)
if n_elements(alpha) eq 0 then alpha = dblarr(nstars)

if n_elements(alpha) ne nstars then message, 'TEFF and ALPHA must have the same number of elements'
if n_elements(logg) ne nstars then message, 'TEFF and LOGG must have the same number of elements'
if n_elements(feh) ne nstars then message, 'TEFF and FEH must have the same number of elements'
if n_elements(av) ne nstars then message, 'TEFF and AV must have the same number of elements'
if n_elements(distance) ne nstars then message, 'TEFF and DISTANCE must have the same number of elements'
if n_elements(lstar) ne nstars then message, 'TEFF and LSTAR must have the same number of elements'
;if n_elements(errscale) ne nstars then message, 'TEFF and ERRSCALE must have the same number of elements'
if n_elements(specphotpath) eq 0 then specphotpath=''

;; store the atmosphere grids in memory to avoid expensive disk IO at each step
common multised_block, wavelength, nwaves, kapv, kapp1, mag, errmag, flux, errflux, blend, sedbands, weff, widtheff, nbands, filter_curves, nspecfiles, labels, spectrophotometry, specblend, filter_curve_sum, zero_point=zero_point

;c=2.9979e14 ;; um/s
if n_elements(pc) eq 0 then pc=3.0857e18 ;; cm
if n_elements(rsun) eq 0 then rsun=6.96e10 ;; cm

;; don't do this every time
if n_elements(flux) eq 0 or keyword_set(redo) then begin
   nwaves = 24000
   wavelength = findgen(nwaves)/1000+0.1 ;; wavelength scale on which to interpolate

   ;; read in the SED bands
   readsedfile, sedfile, nstars, sedbands=sedbands, mag=mag,errmag=errmag,blend=blend,$
                filter_curves=filter_curves, weff=weff, widtheff=widtheff, zero_point=zero_point, $
                flux=flux, errflux=errflux, filter_curve_sum=filter_curve_sum

   blend0 = blend

   nbands=n_elements(weff)
   
   readcol,filepath('extinction_law.ascii',root_dir=getenv('EXOFAST_PATH'),subdir='sed'),klam,kkap,/silent
   kapv = interpol(kkap,klam,0.55)
   kapp1 = interpol(kkap,klam,wavelength)

   ;; if we want to do spectrophotometry, read it in here
   if specphotpath ne '' then begin
      specfiles = file_search(specphotpath,count=nspecfiles)
      specblend = bytarr(nspecfiles,nstars)
      spectrophotometry = ptrarr(nspecfiles,/allocate_heap)
      labels = strarr(nspecfiles)
      for i=0L, nspecfiles-1 do begin
         ;; this file must be um, erg/s/cm^2, erg/s/cm^2
         readcol, specfiles[i], specwave, specflux, specfluxerr,format='d,d,d',/silent
         *spectrophotometry[i] = [[specwave],[specflux],[specfluxerr]]

         fileparts = strsplit(specfiles[i],'.',/extract)
         labels[i] = fileparts[1]
         
         ;; the filename can encode the blending
         if n_elements(fileparts) ge 3 then begin
            blendtxt = fileparts[2]
            if total(valid_num( blendtxt)) eq n_elements(blendtxt) then begin
               blendndx = long(blendtxt)
               if max(blendndx) lt nstars and min(blendndx) ge 0 then begin
                  specblend[blendndx,i] = 1B
               endif
            endif else specblend[*,i] = 1B
         endif else specblend[*,i] = 1B

      endfor
   endif else nspecfiles=0

endif

sed = dblarr(nstars,nwaves)
for j=0L, nstars-1 do begin
   if keyword_set(oned) then begin
      ;; round to 0.5 dex in logg and met rather than interpolate
      logg1 = double(round(logg[j]*2d0))/2d0
      exofast_interp_model,teff[j],logg1[j],feh[j],wavelength,lamflam1temp,alpha=alpha[j],/next, interpfiles=interpfiles
   endif else begin
      exofast_interp_model3d,teff[j],logg[j],feh[j],wavelength,lamflam1temp,alpha=alpha[j],/next, interpfiles=interpfiles, logname=logname, verbose=verbose
   endelse

   ;; interpolation failed, skip
   if ~finite(lamflam1temp[0]) then return, !values.d_infinity

   ;; convert to observed flux
   lamflam1=lamflam1temp*rstar[j]*rstar[j]*rsun*rsun/distance[j]/distance[j]/pc/pc 
   taul1 = kapp1/kapv/1.086*Av[j]
   extinct1 = exp(-taul1)
   sed[j,*] = lamflam1*extinct1
endfor

;; compute the blended flux in each band
modelfluxpos = dblarr(nbands)
modelfluxneg = dblarr(nbands)
if nspecfiles gt 0 then begin
   specphotflux = ptrarr(nspecfiles,/allocate_heap)
endif

for j=0L, nstars-1 do begin
   for i=0L, nbands-1 do begin
      if blend[i,j] eq 1 then modelfluxpos[i] += total(sed[j,*]*filter_curves[i,*])/filter_curve_sum[i] $
      else if blend[i,j] eq -1 then modelfluxneg[i] += total(sed[j,*]*filter_curves[i,*])/filter_curve_sum[i]
   endfor
   for i=0L, nspecfiles-1 do begin
      ;; interpolate model flux onto observed scale
      if abs(specblend[i,j]) eq 1 then intflux = interpol(sed[j,*], wavelength, (*spectrophotometry[i])[*,0])
      if specblend[i,j] eq 1 then begin
         if n_elements((*specphotflux[i])) eq 0 then (*specphotflux[i]) = intflux $
         else (*specphotflux[i]) += intflux
      endif
      if specblend[i,j] eq -1 then begin
         if n_elements((*specphotflux[i])) eq 0 then (*specphotflux[i]) = -intflux $
         else (*specphotflux[i]) -= intflux
      endif
   endfor
endfor

;stop

sedchi2=0d0

;; chi2 from broad band photometry
relative = where(modelfluxneg ne 0,complement=absolute)

;; gross. there has to be a more elegant way...
if absolute[0] ne -1 then sedchi2 += exofast_like(flux[absolute]-modelfluxpos[absolute],0d0,errflux[absolute]*errscale,/chi2)
if relative[0] ne -1 then sedchi2 += exofast_like(mag[relative]+2.5d0*alog10(modelfluxpos[relative]/modelfluxneg[relative]),0d0,errmag[relative]*errscale,/chi2)

;; chi2 from spectrophotometry
for i=0L, nspecfiles-1 do begin
   sedchi2 += exofast_like((*specphotflux[i])-(*spectrophotometry[i])[*,1],0d0,(*spectrophotometry[i])[*,2]*sperrscale[i],/chi2)
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
   
   ;; define the limits; plot the model flux
   xmin = min(weff, max=xmax)
   xmax = 30
   xmin = 0.1
   ymin = alog10(min([flux[absolute],flux[absolute]-errflux[absolute]])) ;,reform(atmospheres,n_elements(atmospheres))]))
   ymax = alog10(max([flux[absolute],flux[absolute]+errflux[absolute],total(sed,1)]))

   if finite(range[0]) then xmin = range[0]
   if finite(range[1]) then xmax = range[1]
   if finite(range[2]) then ymin = range[2]
   if finite(range[3]) then ymax = range[3]
   
   plot, [0], [0], /xlog, ytitle=ytitle, yrange=[ymin,ymax], xrange=[xmin,xmax], /xs, position=position1, xtickformat='(A1)'
      
   ;; set colors for each star and all stars
   pointcolors = dblarr(nbands) + colors[0]

   if nstars gt 1 then begin
      ;; plot atmosphere for each individual star
      for j=0L, nstars-1 do begin
         oplot, wavelength, alog10(smooth(sed[j,*],10)), color=colors[(j+1) mod ncolors]
      endfor
      
      
      for i=0L, nbands-1 do begin      
         if total(blend[i,*] eq 1) eq 1d0 then pointcolors[i] = colors[((where(blend[i,*] eq 1L))[0]+1) mod ncolors] ;$ 
;         else if total(blend[i,*] eq 1) eq nstars then pointcolors[i] = colors[0] $ ;; all 
;         else pointcolors[i] = colors[((where(blend[i,*] eq -1L))[-1]+1) mod ncolors] ;; relative, use color of fainter star
      endfor

      ;; plot atmospheres for each combination of blended stars 
      ;; except 1 and all, handled separately
      starnames = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
      legendlabels = starnames[lindgen(nstars)]
      legendcolors = colors[(lindgen(nstars)+1) mod ncolors]
      colorndx = nstars+1
      
      for i=0L, nbands-1L do begin 
         ;; if the observed band is some combination of more than one but not all stars
         if total(abs(blend[i,*])) gt 1 and total(blend[i,*]) ne nstars then begin
            starstr = strjoin(starnames[where(blend[i,*] eq 1)],'+')
            is_diffmag = ((where(blend[i,*] eq -1))[0] ne -1)
            if is_diffmag then begin
               starstr += '-' + strjoin(starnames[where(blend[i,*] eq -1)],'-')
            endif
            
            ;; plot blended atmospheres for all supplied combinations
            ;; include differential photometry
            ;; if differential photometry supplied as (A+B) - (C+D), 
            ;; plot A+B and C+D
            for k=-1,1,2 do begin 
               blendstarndx = where(blend[i,*] eq k, nblend)
               if nblend gt 1 and nblend ne nstars then begin
                  legendtxt = strjoin(starnames[blendstarndx],'+')
                  legendndx = where(legendlabels eq legendtxt)
                  if legendndx[0] eq -1 then begin ;; if I haven't done this blend yet
                     blended_atmosphere = total(sed[blendstarndx,*],1)
                     color = colors[colorndx mod ncolors]
                     oplot, wavelength, alog10(smooth(blended_atmosphere,10)), color=color
                     legendlabels = [legendlabels,legendtxt]
                     legendcolors = [legendcolors,color]
                     colorndx++
                  endif
               endif
            endfor
            
            ;; now choose the plot point color, unique to each supplied combination
            match = where(legendlabels eq starstr,nmatch) 
            if nmatch eq 0 then begin
               color = colors[colorndx mod ncolors]
               legendlabels = [legendlabels,starstr]
               legendcolors = [legendcolors,color]
               colorndx++
               pointcolors[i] = color
            endif else pointcolors[i] = legendcolors[match[0]]
         endif 
      endfor ;; each observed band

      legendlabels = [legendlabels,strjoin(starnames[0:nstars-1],'+')]
      legendcolors = [legendcolors,colors[0]]
   
   endif

   ;; plot all stars blended together
   oplot, wavelength, alog10(smooth(total(sed,1),10)),color=colors[0]

   ;; overplot the spectrophotometry data
   if nspecfiles gt 0 then specphotcolors = lonarr(nspecfiles)
   for i=0L, nspecfiles-1 do begin
      color = colors[colorndx mod ncolors]
      oplot, (*spectrophotometry[i])[*,0], alog10((*spectrophotometry[i])[*,1]), color=color
      legendlabels = [legendlabels,labels[i]]
      legendcolors = [legendcolors,color]
      specphotcolors[i] = color
      colorndx++
   endfor

   ;; add the spectrophotometry to the legend
   if nstars gt 1 or nspecfiles gt 0 then exofast_legend, legendlabels, textcolors=legendcolors,/top,/right,charsize=0.75

   ;; for the residuals in the lower panel
   residuals = dblarr(nbands)
   res_errhi = dblarr(nbands)
   res_errlo = dblarr(nbands)
   
   for i=0, nbands-1 do begin

      ;; plot model bands (blue filled circles)
      relative = where(blend[i,*] eq -1)
      if relative[0] eq -1 then begin
         oplot, [weff[i]],[alog10(modelfluxpos[i])], psym=8 ;; blue points
         residuals[i] = (flux[i] - (modelfluxpos[i]-modelfluxneg[i]))/errflux[i]
         res_errhi[i] = (flux[i] - (modelfluxpos[i]-modelfluxneg[i])+errflux[i])/errflux[i]
         res_errlo[i] = (flux[i] - (modelfluxpos[i]-modelfluxneg[i])-errflux[i])/errflux[i]
      endif else begin
         ;; this is relative photometry
         ;; only plot the positive stars' model flux
         oplot, [weff[i]],[alog10(modelfluxpos[i])],psym=8 ;; blue points

         ;; pos = neg + deltamag

         ;; overwrite global flux value of relative fluxes with neg star flux + deltamag
         ;; add the model for the positive one 
         ;; to plot relative data on an absolute scale
         flux[i] = modelfluxneg[i]*10^(-0.4*(mag[i])) ;; equal to positive stars' flux
         errflux[i] = flux[i]*alog(10d0)/2.5d0*errmag[i]

         ;; gross. There has to be a more elegant way...
         residuals[i] = (-2.5d0*alog10(modelfluxpos[i]/modelfluxneg[i])-mag[i])/(errmag[i]*errscale)
         res_errhi[i] = (-2.5d0*alog10(modelfluxpos[i]/modelfluxneg[i])-mag[i]+errmag[i]*errscale)/(errmag[i]*errscale)
         res_errlo[i] = (-2.5d0*alog10(modelfluxpos[i]/modelfluxneg[i])-mag[i]-errmag[i]*errscale)/(errmag[i]*errscale)

      endelse

      ;; plot the observed bands (red points with 2D error bars)
      ;; oploterror has too many dependencies; do it myself
      ;; x error bar (red points)
      oplot, [weff[i]-widtheff[i]/2d0,weff[i]+widtheff[i]/2d0], alog10([flux[i],flux[i]]), color=colors[1]
      ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(weff[i]-widtheff[i]/2d0,alog10(flux[i]),/to_device)
      xy2 = convert_coord(weff[i]+widtheff[i]/2d0,alog10(flux[i]),/to_device)
      plots, [xy1[0],xy1[0]], [xy1[1]-ebw,xy1[1]+ebw], color=colors[1],/device
      plots, [xy2[0],xy2[0]], [xy2[1]-ebw,xy2[1]+ebw], color=colors[1],/device
      
      ;; y error bar (red points)
      oplot, [weff[i],weff[i]], alog10([flux[i]-errflux[i],flux[i]+errflux[i]]), color=colors[1]
      ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(weff[i],alog10(flux[i]-errflux[i]),/to_device)
      xy2 = convert_coord(weff[i],alog10(flux[i]+errflux[i]),/to_device)
      plots, [xy1[0]-ebw,xy1[0]+ebw], [xy1[1],xy1[1]], color=colors[1],/device
      plots, [xy2[0]-ebw,xy2[0]+ebw], [xy2[1],xy2[1]], color=colors[1],/device
      
   endfor
   
   ;; now do the residuals (lower panel)
   ;; round to 0.5
   ymin = floor(min(res_errlo)/0.5)*0.5
   ymax = ceil(max(res_errhi)/0.5)*0.5
   
   ;; make yrange symmetric
   if abs(ymin) gt abs(ymax) then ymax =  abs(ymin)
   if abs(ymax) gt abs(ymin) then ymin = -abs(ymax)
   
   ;; if user supplied plotting ranges, use them
   if finite(range[4]) then ymin = range[4]
   if finite(range[5]) then ymax = range[5]
   
   ;; plot the shell of the lower panel
   plot, [0],[0], position=position2, /noerase, $
         xrange=[xmin,xmax], xtitle=xtitle, /xlog, $
         yrange=[ymin,ymax]/0.7, ytitle=ytitle2,$ ;ytitle='O-C',$;ytitle=textoidl('O-C (\sigma)'), $
         /xstyle, /ystyle, yminor=2,yticks=2,ytickv=[ymin,0d0,ymax]
   ;; plot a dashed line at 0
   oplot, [xmin,xmax], [0d0,0d0], linestyle=1,color=colors[1]

   ;; plot the spectrophotometry residuals
   for i=0L, nspecfiles-1 do begin
      oplot, (*spectrophotometry[i])[*,0],$
             ((*specphotflux[i])-(*spectrophotometry[i])[*,1])/(*spectrophotometry[i])[*,2],$
             color=specphotcolors[i]
   endfor

   for i=0L, nbands-1 do begin
      plotsym, 0, symsize, /fill, color=pointcolors[i]
      ;; plot the data points
      oplot, [weff[i]], [residuals[i]], psym=8;, color=pointcolors[i]
      
      ;; x error bar
      oplot, [weff[i]-widtheff[i]/2d0,weff[i]+widtheff[i]/2d0], [residuals[i],residuals[i]], color=pointcolors[i]
      ebw = !d.y_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(weff[i]-widtheff[i]/2d0,residuals[i],/to_device)
      xy2 = convert_coord(weff[i]+widtheff[i]/2d0,residuals[i],/to_device)
      
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
      
      ;; y error bar on residuals
      oplot, [weff[i],weff[i]], [res_errlo[i],res_errhi[i]], color=pointcolors[i]
      ebw = !d.x_vsize/100d0 ;; error bar width = 1% of device size
      xy1 = convert_coord(weff[i],res_errlo[i],/to_device)
      xy2 = convert_coord(weff[i],res_errhi[i],/to_device)
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
      
;      exofast_forprint, sedbands, weff, widtheff, sed, errflux, flux, flux-flux,startxt, textout=residualfilename, $
;                        comment='# Filtername, Center wavelength (um), half bandpass (um), flux, error, flux, residuals (erg/s/cm^2), star indices', $
;                        format='(a20,x,f0.6,x,f0.6,x,e0.6,x,e0.6,x,e0.6,x,e0.6,x,a)'
      
   endif
   set_plot, mydevice
   
end

return, sedchi2

end

