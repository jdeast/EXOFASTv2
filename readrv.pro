function readrv, filename

label = (strsplit(filename,'.',/extract))(1)

line = ''
openr, lun, filename,  /get_lun
readf, lun, line
if strpos(line,'SB2_') ne -1 then begin
   planet = long((strsplit(line,'SB2_',/extract,/regex))[1])
endif else planet = -1L
free_lun, lun

;; Read the RV data file into the RV structure
;; /QUICK isn't necessary for speed, but helps the user find
;; badly formatted files, especially when multi-threading
readcol, filename, bjd, rv, err, format='d,d,d',/silent,comment='#';,/quick

bi = err + !values.d_nan
bierr = err  + !values.d_nan

return, create_struct('bjd',bjd,'rv',rv,'err',err, 'bi',bi,'bierr',bierr,'label',label,'residuals',dblarr(n_elements(bjd)), 'rm',dblarr(n_elements(bjd)), 'planet',planet)

end
