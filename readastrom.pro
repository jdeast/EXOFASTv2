function readastrom, filename

if not file_test(filename) then message, 'Astrometry file (' + filename + ') does not exist'

basename = file_basename(filename)
if n_elements(strsplit(basename,'.',/extract)) lt 4 then message, 'filename (' + basename + ') must have format EPOCH.BAND.LOCATION.whateveryouwant (see readastrom.pro for details)'

epoch = (strsplit(basename,'.',/extract))[0]
location = (strsplit(basename,'.',/extract))[2]

;; Read the transit data file into a structure
;; (with an arbitary number of detrending variables)
band = (strsplit(basename,'.',/extract))(1)
if band eq 'Sloanu' then begin
   band = 'Sloanu'
   bandname = "u'"
endif else if band eq 'Sloang' then begin
   band = 'Sloang'
   bandname = "g'"
endif else if band eq 'Sloanr' then begin
   band = 'Sloanr'
   bandname = "r'"
endif else if band eq 'Sloani' then begin
   band = 'Sloani'
   bandname = "i'"
endif else if band eq 'Sloanz' then begin
   band = 'Sloanz'
   bandname = "z'"
endif else bandname = band

allowedbands = ['U','B','V','R','I','J','H','K',$
                'Sloanu','Sloang','Sloanr','Sloani','Sloanz',$
                'Kepler','TESS','CoRoT','Spit36','Spit45','Spit58','Spit80',$
                'u','b','v','y']

if (where(allowedbands eq band))[0] eq -1 then message, 'Filter (' + band + ') not allowed'



line = ''
openr, lun, filename, /get_lun
readf, lun, line
entries = strsplit(line, /extract)
free_lun, lun
nentries = n_elements(entries)
if nentries eq 5 then begin 
   if location eq 'Earth' or location eq 'Gaia' or location eq 'Hipparcos' then $
      readcol, filename, bjdtdb, ra, dec, raerr, decerr, format='d,d,d,d,d',/silent $
   else message, 'Location not recognized. Must supply observatory position in file'
   
   if location eq 'Earth' then begin
      tbase = 0d0 ;; do I want to allow the user to specify tbase?
      ephemfile = find_with_def('JPLEPH.405','ASTRO_DATA')
      mindate = min(bjdtdb+tbase,max=maxdate) - 1.d0 & maxdate+=1
      ;; use the geocenter of the Earth
      JPLEPHREAD,ephemfile,pinfo,pdata,[mindate,maxdate], $
                 status=status, errmsg=errmsg
      JPLEPHINTERP,pinfo,pdata,bjdtdb,x_earth,y_earth,z_earth,$
                   /earth,posunits="AU", tbase=tbase     
      obspos = transpose([[x_earth],[y_earth],[z_earth]])

      ;; for plotting a pretty model
      npretty = (maxdate-mindate)*4 ;; every 6 hours
      prettytime = mindate + (maxdate-mindate)*dindgen(npretty)/(npretty-1)
      JPLEPHINTERP,pinfo,pdata,prettytime,x_earth,y_earth,z_earth,$
                   /earth,posunits="AU", tbase=tbase     
      prettyobspos = transpose([[x_earth],[y_earth],[z_earth]])

   endif else if location eq 'Gaia' then begin
      message, 'Gaia ephemeris not yet supported!!'
   endif else if location eq 'Hipparcos' then begin
      message, 'Hippoarcos ephemeris not yet supported!!'
   endif
endif else if nentries eq 8 then begin
   readcol, filename, bjdtdb, ra, dec, raerr, decerr, xobs, yobs, zobs, format='d,d,d,d,d,d,d,d',/silent
   obspos = transpose([[xobs],[yobs],[zobs]])
endif else begin
   message, 'Astrometry file must have either 5 or 8 columns (see readastrom.pro for details)'
endelse
   
astrom=create_struct('bjdtdb',bjdtdb,'radec',transpose([[ra],[dec]]), 'err',transpose([[raerr],[decerr]]/3600d3),$
                     'obspos',obspos,'epoch',epoch,'band',band,'prettytime',prettytime,'prettyobspos',prettyobspos,$
                     'label', location, 'residuals',transpose([[ra],[dec]])*0d0,'model',transpose([[ra],[dec]])*0d0)

return, astrom

end
