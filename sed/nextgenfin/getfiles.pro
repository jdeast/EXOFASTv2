pro getfiles



for i=1L, 13631L do begin

   spawn, 'curl "http://svo2.cab.inta-csic.es/theory/newov/ssap.php?model=bt-nextgen-agss2009&fid=' + strtrim(i,2) + '&format=ascii"', output
   teff = (strsplit(output[1],/extract))[3]
   logg = (strsplit(output[2],/extract))[3]
   meta = (strsplit(output[3],/extract))[3]
   alpha = (strsplit(output[4],/extract))[3]

   npoints = n_elements(output)-10
   wavelength = dblarr(npoints)
   flux = dblarr(npoints)

   for j=0L, npoints-1 do begin
      arr = strsplit(output[j+9],/extract)
      wavelength[j] = arr[0]
      flux[j] = arr[1]
   endfor

   outname = string(round(double(teff)/100.0), double(logg), double(meta), format='("lte",i2,"-",f3.1,"-",f3.1,".NextGen.spec")')

   openw, lun, outname, /get_lun
   printf, lun, teff, logg, meta, alpha 'Teff, logg, [M/H], alpha'
   printf, lun, npoints, 'number of wavelength points'
   printf, lun, wavelength
   printf, flux
   free_lun, lun

   stop

endfor


;readcol, 'filelist.txt', files, format='a'
;for i=0L, nfiles-1 do begin
;
;   extractedfile = (strsplit(files[i],'.gz',/extract,/regex))[0]
;   if ~file_test(files[i]) and ~file_test(extractedfile) then begin
;      arr = strsplit(files[i],'-.',/extract)
;      zfolder = 'Z-' + arr[3] + '.' + arr[4] 
;
;      spawn, 'wget ftp://calvin.physast.uga.edu/pub/NextGen/Spectra/' + zfolder + '/' files[i]
;      
;   endif
;
;endfor

end
