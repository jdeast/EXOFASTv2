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


;; BLEND is an NSTARSxNBANDS boolean specifying which stars each magnitude applies to

function mistmultised, teff, logg, feh, av, distance, lstar, errscale, sedfile, redo=redo, psname=psname, debug=debug, atmospheres=atmospheres, wavelength=wavelength, logname=logname,range=range,blend0=blend0

nstars = n_elements(teff)
if n_elements(logg) ne nstars then message, 'TEFF and LOGG must have the same number of elements'
if n_elements(feh) ne nstars then message, 'TEFF and FEH must have the same number of elements'
if n_elements(av) ne nstars then message, 'TEFF and AV must have the same number of elements'
if n_elements(distance) ne nstars then message, 'TEFF and DISTANCE must have the same number of elements'
if n_elements(lstar) ne nstars then message, 'TEFF and LSTAR must have the same number of elements'
if n_elements(errscale) ne nstars then message, 'TEFF and ERRSCALE must have the same number of elements'

;; store the BC tables in memory to avoid expensive disk IO at each step
common BC_block, bcarrays, teffgrid, logggrid, fehgrid, avgrid, sedbands, mags, errs, filterprops, blend

if n_elements(teffgrid) eq 0 or keyword_set(redo) then begin

   ;; read in the file with observed magnitudes
   line = ''
   nlines = file_lines(sedfile)
   sedbands = strarr(nlines)
   mags = dblarr(nlines)
   errs = dblarr(nlines)+99d0
   blend = intarr(nlines,nstars)

   openr, lun, sedfile, /get_lun
   for i=0L, nlines-1 do begin
      readf, lun, line
      line = (strsplit(line,'#',/extract,/preserve_null))[0] ;; remove comments
      ;if strpos(line,'#') eq 0 then continue
      entries = strsplit(line,/extract)
      if n_elements(entries) lt 3 then begin
         continue ;; not a legal line
      endif else begin
         sedbands[i] = entries[0]
         mags[i] = double(entries[1])
         errs[i] = double(entries[2])
      endelse
      if n_elements(entries) eq 5 then begin
         if strpos(entries[4],'-') ne -1 then begin
            ;; handle differential magnitudes
            posndx = long(strsplit((strsplit(entries[4],'-',/extract))[0],',',/extract))

            good = where(posndx lt nstars and posndx ge 0,complement=bad)
            if bad[0] ne -1 then printandlog, 'WARNING: STARNDX (' + strtrim(posndx[bad],2) + ') in SEDFILE (' + sedfile +') does not correspond to a star, ignoring', logname
            if good[0] eq -1 then begin
               printandlog, 'WARNING: No good positive STARNDX values in SEDFILE (' + sedfile +') line: ' + line, logname
               continue
            endif
            blend[i,posndx[good]] = 1

            negndx = long(strsplit((strsplit(entries[4],'-',/extract))[1],',',/extract))
            good = where(negndx lt nstars and negndx ge 0,complement=bad)
            if bad[0] ne -1 then printandlog, 'WARNING: STARNDX (' + strtrim(negndx[bad],2) + ') in SEDFILE (' + sedfile +') does not correspond to a star, ignoring', logname
            if good[0] eq -1 then begin
               printandlog, 'WARNING: No good negative STARNDX values in SEDFILE (' + sedfile +') line: ' + line, logname
               continue
            endif
            blend[i,negndx[good]] = -1

         endif else begin
            ;; handle unblended magnitudes
            starndx = long(strsplit(entries[4],',',/extract))
            good = where(starndx lt nstars and starndx ge 0,complement=bad)
            if bad[0] ne -1 then printandlog, 'WARNING: STARNDX (' + strtrim(starndx[bad],2) + ') in MISTSEDFILE (' + sedfile +') does not correspond to a star, ignoring', logname
            if good[0] eq -1 then begin
               printandlog, 'WARNING: No good STARNDX values in MISTSEDFILE (' + sedfile +') line: ' + line, logname
               continue
            endif
            blend[i,starndx[good]] = 1
         endelse
      endif else begin
         ;; by default, assume all stars are blended together
         blend[i,*] = 1
      endelse
   endfor
   free_lun, lun

;   readcol, sedfile, sedbands, mags, errs, format='a,d,d', comment='#', /silent
   good = where(errs lt 1d0, ngood)
   if ngood gt 1 then begin
      sedbands = sedbands[good]
      mags = mags[good]
      errs = errs[good]
      blend = blend[good,*]
  endif else begin
      print, 'Bands must have errors less than 1 mag; no good bands'
      stop
   endelse
   nbands = n_elements(sedbands)

;   readcol, filepath('filternames.txt', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist']), keivanname, mistname, format='a,a', comment='#',/silent
   readcol, filepath('filternames2.txt', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist']), keivanname, mistname, claretname, svoname, format='a,a,a,a', comment='#',/silent

   ;; this is the dimension for each axis of the array
   restore, filepath('mist.sed.grid.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
   for i=0L, nbands-1 do begin

      filename = filepath(sedbands[i] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
      if not file_test(filename) then begin
         match = (where(keivanname eq sedbands[i]))[0]
         filename = filepath(mistname[match] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
         if match[0] eq -1 or not file_test(filename) then begin
            match = (where(svoname eq sedbands[i]))[0]
            filename = filepath(mistname[match] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
            if match[0] eq -1 or not file_test(filename) then begin
               print, sedbands[i] + ' not supported, remove or comment out from ' + sedfile
               stop
            endif
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

nbands = n_elements(sedbands)
bcs=dblarr(nbands,nstars)

for j=0L, nstars-1 do begin
   
   ;; bolometric corrections
   teff_ndx = getgridpoint(teffgrid,teff[j])
   logg_ndx = getgridpoint(logggrid,logg[j])
   feh_ndx = getgridpoint(fehgrid,feh[j])
   av_ndx = getgridpoint(avgrid,av[j])

   ;; quadrilinear interpolation of the MIST Bolometric Correction tables (Teff, logg, [Fe/H], A_V)  
   for i=0L, nbands-1 do $ 
      bcs[i,j] = ninterpolate(bcarrays[*,*,*,*,i],[teff_ndx,logg_ndx,feh_ndx,av_ndx])

endfor

blend0 = blend

modelmag = -2.5d0*alog10(lstar##replicate(1d0,nbands)) + 4.74d0 - bcs + mu##replicate(1d0,nbands)
modelflux = 10^(-0.4*modelmag)
magresiduals = dblarr(nbands)
if nstars eq 1 then begin
   blendmag = modelmag
   blendflux = modelflux
   magresiduals = mags-modelmag
endif else begin 
   blendfluxpos = total(modelflux*(blend gt 0), 2)
   blendfluxneg = total(modelflux*(blend lt 0), 2)
   blendfluxneg[where(blendfluxneg eq 0d0)] = 1d0
   blendmag = -2.5d0*alog10(blendfluxpos/blendfluxneg)
   magresiduals = mags-blendmag
   blendflux = total(modelflux*blend,2)
endelse

;; calculate the likelihood
;; sedchi2 = exofast_like(mags-blendmag,0d0,errs*errscale[0],/chi2)
sedchi2 = exofast_like(magresiduals,0d0,errs*errscale[0],/chi2)

if 0 then begin
forprint, modelflux[*,0], modelflux[*,1], mags,blendmag,((mags-blendmag)/(errs*errscale[0]))^2, /textout
print, 'lstar = ' + strtrim(lstar,2)
print, 'teff = ' + strtrim(teff,2)
print, 'logg = ' + strtrim(logg,2)
print, 'feh = ' + strtrim(feh,2)
print, 'av = ' + strtrim(av,2)
print, 'chi2 = ' + strtrim(sedchi2,2)
;wait, 0.1
;stop
endif

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
  
   wp = transpose(filterprops.lambda_eff)/1d4 ;; effective wavelength, in um
   widthhm = transpose(filterprops.w_eff)/1d4/2d0 ;; half width, in um
   
   ;; convert everything to AB mag scale
   vega = where(filterprops.type eq 'Vega')
   abmags = mags
   modelabmags = blendmag;modelmags
   if vega[0] ne -1 then begin
      abmags[vega] += filterprops[vega].vegaab
      modelabmags[vega] += filterprops[vega].vegaab
   endif

   ;; compute the absolute fluxes
   zp = 3631d0*3d-9/wp
   flux = zp*10^(-0.4d0*abmags)
   modelblendflux = zp*10^(-0.4d0*modelabmags)
   ;; correct flux array for relative fluxes
   for i=0L, nbands-1 do begin
      relative = where(blend[i,*] eq -1)
      if relative[0] ne -1 then begin
         ;; assume model flux for positive stars
         blendmagpos = -2.5d0*alog10(blendfluxpos[i])
         blendmagneg = -2.5d0*alog10(blendfluxneg[i])
         if vega[i] ne -1 then begin
            blendmagpos += filterprops[i].vegaab
            blendmagneg += filterprops[i].vegaab
         endif
         flux[i] = zp[i]*10^(-0.4d0*(blendmagpos-mags[i]))
         modelblendflux[i] = zp[i]*10^(-0.4d0*blendmagneg) ;; only plot the fainter star
      endif
   endfor
   fluxerr = flux*alog(10d0)/2.5d0*errs

   ;; define the limits; plot the model flux
   xmin = min(wp, max=xmax)
   xmax = 30
   xmin = 0.1
   ymin = alog10(min([modelblendflux,flux-fluxerr])) ;,reform(atmospheres,n_elements(atmospheres))]))
   ymax = alog10(max([modelblendflux,flux+fluxerr,total(atmospheres,1)]))
   if finite(range[0]) then xmin = range[0]
   if finite(range[1]) then xmax = range[1]
   if finite(range[2]) then ymin = range[2]
   if finite(range[3]) then ymax = range[3]

   plot, [0], [0], /xlog, ytitle=ytitle, yrange=[ymin,ymax], xrange=[xmin,xmax], /xs, position=position1, xtickformat='(A1)'

   ;; plot model atmospheres for each individual star
   for i=0L, nstars-1 do begin
      oplot, wavelength, alog10(smooth(atmospheres[i,*],10)), color=colors[(i+1) mod ncolors]
   endfor

   ;; set colors for each star and all stars
   pointcolors = dblarr(nbands)
   for i=0L, nbands-1 do begin
      if total(blend[i,*] eq 1) eq 1d0 then pointcolors[i] = colors[((where(blend[i,*] eq 1L))[0]+1) mod ncolors] $ ;; only one star, use its color
      else if total(blend[i,*] eq 1) eq nstars then pointcolors[i] = colors[0]; $ ;; all blended, use black
;      else pointcolors[i] = colors[((where(blend[i,*] eq -1L))[-1]+1) mod ncolors] ;; relative, use new colo
   endfor

   ;; plot model atmospheres for each combination of blended stars 
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
                  blended_atmosphere = total(atmospheres[blendstarndx,*],1)
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


if 0 then begin
   for i=0L, nbands-1L do begin
      if total(abs(blend[i,*])) gt 1 and total(blend[i,*]) ne nstars then begin
         ; for j=0L, i-1 do begin
         ;; if no other band before has this combination of stars
         starstr = strjoin(starnames[where(blend[i,*] eq 1)],'+')
         is_diffmag = ((where(blend[i,*] eq -1))[0] ne -1)
         if is_diffmag then begin
            starstr += '-' + strjoin(starnames[where(blend[i,*] eq -1)],'-')
         endif
         match = where(legendlabels eq starstr,nmatch)
         if nmatch eq 0 then begin
            ;; plot blended atmosphere if not diffmag
            if ~is_diffmag then begin
               blended_atmosphere = transpose(atmospheres[0,*])*0d0
               for k=0L, nstars-1 do begin
                  if blend[i,k] then blended_atmosphere += atmospheres[k,*]
               endfor
               oplot, wavelength, alog10(smooth(blended_atmosphere,10)), color=colors[colorndx mod ncolors]
            endif
            color = colors[colorndx mod ncolors]
            legendlabels = [legendlabels,starstr]
            legendcolors = [legendcolors,color]
            colorndx++
            pointcolors[i] = color
         endif else pointcolors[i] = legendcolors[match[0]]
         ; endfor
      endif
   endfor
endif

   ;; plot all model atmospheres blended together
   oplot, wavelength, alog10(smooth(total(atmospheres,1),10)),color=colors[0]
   legendlabels = [legendlabels,strjoin(starnames[0:nstars-1],'+')]
   legendcolors = [legendcolors,colors[0]]
   if nstars gt 1 then exofast_legend, legendlabels, textcolors=legendcolors,/top,/right,charsize=0.75

   ;; plot bands
   oplot, wp, alog10(modelblendflux), psym=8

   ;; oploterror has too many dependencies; do it myself
   for i=0, nbands-1 do begin
      ;; plot model bands (blue filled circles)
      relative = where(blend[i,*] eq -1)
      if relative[0] eq -1 then begin
         ;; plot all stars
         oplot, [wp[i]], [alog10(modelflux[i])], psym=8
      endif else begin
         ;; for relative photometry, only plot the fainter star's model
         oplot, [wp[i,relative]], [alog10(modelflux[i,relative])], psym=8
      endelse

      ;; plot observed bands (red 2d error bars)
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
   residuals = (flux-modelblendflux)/fluxerr
   res_errhi = (flux-modelblendflux+fluxerr)/fluxerr
   res_errlo = (flux-modelblendflux-fluxerr)/fluxerr

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
         yrange=[ymin,ymax]/0.7, ytitle=ytitle2,$;ytitle='O-C',$;ytitle=textoidl('O-C (\sigma)'), $
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

      exofast_forprint, filterprops.name, wp, widthhm, flux, fluxerr, modelblendflux, flux-modelblendflux,startxt, textout=residualfilename, $
                        comment='# Filtername, Center wavelength (um), half bandpass (um), flux, error, modelflux, residuals (erg/s/cm^2), star indices', $
                        format='(a20,x,f0.6,x,f0.6,x,e0.6,x,e0.6,x,e0.6,x,e0.6,x,a)'

   endif
   set_plot, mydevice
endif


return, sedchi2

end
