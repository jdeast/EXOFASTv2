function readrv, filename

label = (strsplit(filename,'.',/extract))(1)

;; Read the RV data file into the RV structure
readcol, filename, bjd, rv, err, format='d,d,d',/silent,comment='#'

bi = err + !values.d_nan
bierr = err  + !values.d_nan

return, create_struct('bjd',bjd,'rv',rv,'err',err, 'bi',bi,'bierr',bierr,'label',label,'residuals',dblarr(n_elements(bjd)), 'rm',dblarr(n_elements(bjd)))

end
