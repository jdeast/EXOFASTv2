;; linearly interpolate the grid to the value
function getgridpoint, grid, value, box=box

  ngrid = n_elements(grid)

  ;; find the index of the nearest point
  junk = min(abs(grid - value),match)

  if match eq ngrid-1 then begin
     ;; pegged up against the end of the grid
     ;; linearly extrapolate from the last two points
     ndx = match + (value-grid[ngrid-1])/(grid[ngrid-1]-grid[ngrid-2])
  endif else if match eq 0 then begin
     ;; pegged up against the beginning of the grid
     ;; linearly extrapolate from the first two points
     ndx = match + (value-grid[0])/(grid[1]-grid[0])
  endif else begin
     ;; in the middle of the grid
     ;; linearly interpolate from the two nearest points
     if value gt grid[match] then begin
        ndx = match + (value-grid[match])/(grid[match+1]-grid[match])
     endif else begin
        ndx = match + (value-grid[match])/(grid[match]-grid[match-1])
     endelse
  endelse

  return, ndx

end

function mistsed, teff, logg, feh, av, distance, lstar, errscale, sedfile, redo=redo, psname=psname, debug=debug

;; store the BC tables in memory to avoid expensive disk IO at each step
common BC_block, bcarrays, teffgrid, logggrid, fehgrid, avgrid, bands, mags, errs, filterprops

if n_elements(teffgrid) eq 0 or keyword_set(redo) then begin

   ;; read in the file with observed magnitudes
   readcol, sedfile, bands, mags, errs, format='a,d,d', comment='#', /silent
   good = where(errs lt 1d0, ngood)
   if ngood gt 1 then begin
      bands = bands[good]
      mags = mags[good]
      errs = errs[good]
   endif else begin
      print, 'Bands must have errors less than 1 mag; no good bands'
      stop
   endelse
   nbands = n_elements(bands)

   readcol, filepath('filternames.txt', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist']), keivanname, mistname, format='a,a', comment='#',/silent

   ;; this is the dimension for each axis of the array
   restore, filepath('mist.sed.grid.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
   for i=0L, nbands-1 do begin

      filename = filepath(bands[i] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
      if not file_test(filename) then begin
         match = (where(keivanname eq bands[i]))[0]
         filename = filepath(mistname[match] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
         if match[0] eq -1 or not file_test(filename) then begin
            print, bands[i] + ' not supported, remove or comment out from ' + sedfile
            stop
         endif        
      endif  
      restore, filename
      if i eq 0L then begin
         sz = size(bcarray)
         bcarrays = dblarr(sz[1],sz[2],sz[3],sz[4],nbands) ;; bolometric correction arrays
         filterprops = filterproperties ;; meta data for each filter (effective wavelength, effective width, zero point, magnitude system)
      endif else filterprops = [filterprops,filterproperties]
      bcarrays[*,*,*,*,i] = bcarray
   endfor
endif

;; distance modulus
mu = 5d0*alog10(distance)-5d0                   

;; bolometric corrections
teff_ndx = getgridpoint(teffgrid,teff)
logg_ndx = getgridpoint(logggrid,logg)
feh_ndx = getgridpoint(fehgrid,feh)
av_ndx = getgridpoint(avgrid,av)

;; quadrilinear interpolation of the MIST Bolometric Correction tables (Teff, logg, [Fe/H], A_V)  
nbands = n_elements(bands)
bcs=dblarr(nbands)
for i=0L, nbands-1 do $
   bcs[i] = ninterpolate(bcarrays[*,*,*,*,i],[teff_ndx,logg_ndx,feh_ndx,av_ndx])

;; convert BC to observed magnitude
modelmags = -2.5d0*alog10(lstar) + 4.74d0 - bcs + mu

;; calculate the likelihood
sedchi2 = exofast_like(mags-modelmags,0d0,errs*errscale,/chi2)

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

   wp = transpose(filterprops.lambda_eff)/1d4 ;; effective wavelength, in um
   widthhm = transpose(filterprops.w_eff)/1d4/2d0 ;; half width, in um
   
   ;; convert everything to AB mag scale
   vega = where(filterprops.type eq 'Vega')
   abmags = mags
   modelabmags = modelmags
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
   ymin = alog10(min([modelflux,flux-fluxerr]))
   ymax = alog10(max([modelflux,flux+fluxerr]))

   plot, [0], [0], /xlog, xtitle=xtitle,ytitle=ytitle, yrange=[ymin,ymax], xrange=[xmin,xmax], /xs;,/ys
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

   ;; clean up the postscript device
   if keyword_set(psname) then begin
      !p.font=-1
      !p.multi=0
      device, /close
      device, encapsulated=0

      ;; create a residual file
      residualfilename = file_dirname(psname) + path_sep() + file_basename(psname,'.eps') + '.residuals.txt'
      exofast_forprint, filterprops.name, wp, modelflux, flux, fluxerr, fluxerr, flux-modelflux, textout=residualfilename, comment='# Filtername, Wavelength (um), model flux, flux, error, residuals (erg/s/cm^2)' 

   endif
   set_plot, mydevice
endif

return, sedchi2

end
