;; this takes the raw model flux files (retrieved with getfiles.pro)
;; into the required format for exofast_interp_model.pro. Restricting
;; these files to 0.1 to 24 microns reduces the file sizes, load
;; times, and RAM footprint dramatically. Interpolating now saves us
;; time later but may lose generality
pro convertsed

files = file_search(['/home/jeastman/nextgenfin/nextgen.?????.spec',$
                     '/h/onion0/idl/EXOFASTv2/sed/nextgenfin/nextgenfin2/nextgen.?????.spec'],$
                    count=nfiles)

w1 = findgen(24000)/1000+0.1 ;; wavelength scale on which to interpolate (microns)

line = ''
for i=0L, nfiles-1 do begin
   
   openr, lun, files[i], /get_lun
   readf, lun, line ;; discard uselss header
   readf, lun, line
   teff = double((strsplit(line,/extract))[3])
   readf, lun, line
   logg = double((strsplit(line,/extract))[3])
   readf, lun, line
   meta = double((strsplit(line,/extract))[3])
   readf, lun, line
   alpha = double((strsplit(line,/extract))[3])
   free_lun,lun

   teffstr = strtrim(round(teff/100.0),2)
   
   if logg lt 0 then loggsign = '-'$
   else loggsign = '+'
   
   if meta lt 0 then zsign = '-'$
   else zsign = '+'
   
   if alpha lt 0 then asign = '-'$
   else asign = '+'
   
   fmt = '("lte",a,a,f3.1,a,f3.1,a,f3.1,".NextGen.spec.idl")' 
   outname = string(teffstr, loggsign, abs(logg), zsign, abs(meta), asign, abs(alpha), format=fmt)

   if not file_test(outname) then begin
      readcol, files[i], wavelength, flux, comment='#', format='f,f',/silent
      
      lamflam=flux*wavelength
      lamflam1=interpol(lamflam,wavelength,w1*1d4) ;; erg/s/cm^2

      save, lamflam1, filename=outname

;      openw, lun, outname, /get_lun
;      printf, lun, string(teff,logg,meta,alpha,format='(i,x,f4.1,x,f4.1,x,f4.1," Teff, logg, [M/H], alpha")')
;      printf, lun, strtrim(24000,2) + ' number of wavelength points. lambda = findgen(24000)/1000+0.1 ;; microns'
;      printf, lun, lamflam1
;      free_lun, lun

   endif

   print, 'done with ' + outname + ' (' + files[i] + ')'

endfor

stop

end
