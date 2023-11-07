pro readsedfile, sedfile, nstars, sedbands=sedbands, mag=mag, errmag=errmag, blend=blend, $
                 flux=flux, errflux=errflux, filter_curves=filter_curves, weff=weff, widtheff=widtheff, zero_point=zero_point, download_new=download_new, filter_curve_sum=filter_curve_sum


line = ''
nlines = file_lines(sedfile)
sedbands = strarr(nlines)
mag = dblarr(nlines)
errmag = dblarr(nlines)+99d0
flux = dblarr(nlines)
errflux = dblarr(nlines)+99d0
blend = lonarr(nlines,nstars)
filter_curves = dblarr(nlines,24000)
weff = dblarr(nlines)
widtheff = dblarr(nlines)
zero_point = dblarr(nlines)
filter_curve_sum = dblarr(nlines)

openr, lun, sedfile, /get_lun
for i=0L, nlines-1 do begin
   readf, lun, line
   line = (strsplit(line,'#',/extract,/preserve_null))[0] ;; remove comments
   entries = strsplit(line,/extract)
   if n_elements(entries) lt 3 then begin
      continue ;; not a legal line
   endif else begin
      sedbands[i] = entries[0]
      mag[i] = double(entries[1])
      errmag[i] = double(entries[2])
   endelse

   ;; load the filter properties
   idlfile = filepath(sedbands[i]+'.idl',root_dir=getenv('EXOFAST_PATH'),subdir=['sed','filtercurves']) 

   ;; if not recognized, download it
   if not file_test(idlfile) and keyword_set(download_new) then getfilter, sedbands[i]

   ;; if still not recognized, skip it
   if not file_test(idlfile) then begin
      print, sedbands[i] + ' not recognized; skipping'
      continue
   endif

   restore, idlfile
   filter_curves[i,*] = filter.transmission
   weff[i] = filter.weff
   widtheff[i] = filter.widtheff
   zero_point[i] = filter.zero_point
   filter_curve_sum[i] = total(filter_curves[i,*])

   flux[i] = zero_point[i]*10^(-0.4d0*mag[i]);*weff[i]/filter.wavelength
   errflux[i] = flux[i]*alog(10d0)/2.5d0*errmag[i]

   if n_elements(entries) eq 5 then begin

      ;; we can handle differential fluxes
      if strpos(entries[4],'-') ne -1 then begin
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
         starndx = long(strsplit(entries[4],',',/extract))
         good = where(starndx lt nstars and starndx ge 0,complement=bad)
         if bad[0] ne -1 then printandlog, 'WARNING: STARNDX (' + strtrim(starndx[bad],2) + ') in SEDFILE (' + sedfile +') does not correspond to a star, ignoring', logname
         if good[0] eq -1 then begin
            printandlog, 'WARNING: No good STARNDX values in SEDFILE (' + sedfile +') line: ' + line, logname
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

good = where(errmag lt 1d0, ngood)
if ngood gt 1 then begin
   sedbands = sedbands[good]
   mag = mag[good]
   errmag = errmag[good]
   flux = flux[good]
   errflux = errflux[good]
   blend = blend[good,*]

   filter_curves = filter_curves[good,*]
   weff = weff[good]
   widtheff = widtheff[good]
   zero_point = zero_point[good]
   filter_curve_sum = filter_curve_sum[good]

endif else begin
   print, 'Bands must have errors less than 1 mag; no good bands'
   stop
endelse

end


