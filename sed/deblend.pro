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

function deblend, teff, logg, feh, av, distance, lstar, bands, redo=redo, logname=logname

nstars = n_elements(teff)
nbands = n_elements(bands)
if n_elements(logg) ne nstars then message, 'TEFF and LOGG must have the same number of elements'
if n_elements(feh) ne nstars then message, 'TEFF and FEH must have the same number of elements'
if n_elements(av) ne nstars then message, 'TEFF and AV must have the same number of elements'
if n_elements(distance) ne nstars then message, 'TEFF and DISTANCE must have the same number of elements'
if n_elements(lstar) ne nstars then message, 'TEFF and LSTAR must have the same number of elements'

;; store the BC tables in memory to avoid expensive disk IO at each step
common BC_block2, bcarrays, teffgrid, logggrid, fehgrid, avgrid, filterprops

if n_elements(teffgrid) eq 0 or keyword_set(redo) then begin
 
   readcol, filepath('filternames.txt', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist']), keivanname, mistname, claretname, format='a,a,a', comment='#',/silent
   
   ;; this is the dimension for each axis of the array
   restore, filepath('mist.sed.grid.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
   for i=0L, nbands-1 do begin
      

      match = (where(claretname eq bands[i]))[0]
      filename = filepath(mistname[match] + '.idl', root_dir=getenv('EXOFAST_PATH'),subdir=['sed','mist'])
      if match[0] eq -1 or not file_test(filename) then begin
         print, bands[i] + ' not supported, remove or comment out from ' + sedfile
         stop
      endif
      restore, filename ;; this contains bcarray and filterproperties
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
bcs = dblarr(nbands,nstars)

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

modelmag = -2.5d0*alog10(lstar##replicate(1d0,nbands)) + 4.74d0 - bcs + mu##replicate(1d0,nbands)
modelflux = 10^(-0.4*modelmag)

return, modelflux

end
