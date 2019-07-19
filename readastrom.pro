function readastrom, filename

if not file_test(filename) then message, 'Astrometry file (' + filename + ') does not exist'

basename = file_basename(filename)
if n_elements(strsplit(basename,'.',/extract)) lt 4 then message, 'filename (' + basename + ') must have format EPOCH.BAND.LOCATION.whateveryouwant (see readastrom.pro for details)'

epoch = (strsplit(basename,'.',/extract))[0]
if epoch eq 'PA' then begin
   readcol, filename, bjdtdb, rho, rhoerr, pa, paerr, format='d,d,d,d,d',/silent
   userhopa = 1
   ra = 0d0
   dec = 0d0
   raerr = 0d0
   decerr = 0d0
endif else begin
   userhopa = 0
   rho = 0d0
   pa = 0d0
   rhoerr = 0d0
   paerr = 0d0
endelse

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
  
   tbase = 0d0 ;; do I want to allow the user to specify tbase?
   mindate = min(bjdtdb+tbase,max=maxdate) - 1.d0 & maxdate+=1
   if location eq 'Earth' then begin
      ephemfile = find_with_def('JPLEPH.405','ASTRO_DATA')
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

      ;; this requires horizons.exp (distributed with EXOFASTv2) and
      ;; expect (and tcl) to be in your unix path
      jdtdb = [mindate,maxdate]
      baryfile = 'Gaia.' + strtrim(systime(/seconds),2) + '.eph'
      get_eph, jdtdb, 'Gaia', outfile=baryfile

      ;; read in the CSV HORIZONS emphemeris file
      readcol, baryfile, jdlist, junk, xlist, ylist, zlist, $
               format='a,a,d,d,d',delimiter=',',/silent
        
      ;; don't lose precision
      base = double(strmid(jdlist,0,7)) - tbase
      decimal = double(strmid(jdlist,7))
      jdlist = base + decimal
        
      ;; make sure the ephemeris covers the time range
      if min(jdlist) gt min(jd_tdb) or max(jdlist) lt max(jd_tdb) then $
         message,'ERROR: ephemeris does not cover time range -- check ' $
                 + baryfile
      
      ;; interpolate the target positions at each JD
      x_obs = interpol(xlist, jdlist, jd_tdb, /quadratic)
      y_obs = interpol(ylist, jdlist, jd_tdb, /quadratic)
      z_obs = interpol(zlist, jdlist, jd_tdb, /quadratic)
      
      obspos = transpose([[x_obs,y_obs,z_obs]])
      prettyobspos = transpose([[x_list,y_list,z_list]])

      message, 'Gaia ephemeris not yet supported!!'
   endif else if location eq 'Hipparcos' then begin
      message, 'Hipparcos ephemeris not yet supported!!'
   endif
endif else if nentries eq 8 then begin
   readcol, filename, bjdtdb, ra, dec, raerr, decerr, xobs, yobs, zobs, format='d,d,d,d,d,d,d,d',/silent
   obspos = transpose([[xobs],[yobs],[zobs]])
endif else begin
   message, 'Astrometry file must have either 5 or 8 columns (see readastrom.pro for details)'
endelse
   
astrom=create_struct('bjdtdb',bjdtdb,'radec',transpose([[ra],[dec]]), 'err',transpose([[raerr],[decerr]]/3600d3),$
                     'obspos',obspos,'epoch',epoch,'band',band,'prettytime',prettytime,'prettyobspos',prettyobspos,$
                     'label', location, 'residuals',transpose([[ra],[dec]])*0d0,'model',transpose([[ra],[dec]])*0d0, $
                     'userhopa',keyword_set(userhopa), 'rhopa', transpose([[rho],[pa*!dpi/180d0]]), 'rhopaerr', transpose([[rhoerr],[paerr]]))

return, astrom

end
