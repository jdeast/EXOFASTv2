;; read the coefficients from Lindegren+ 2020

pro read_coeffs, filename, delimiter=delimiter, j=j, k=k, g=g, q_jk=q_jk, n=n, m=m

  openr, lun, filename,/get_lun
  line = ''
  readf, lun, line
  j = long(strsplit(line,delimiter,/extract))
  j = j[1:n_elements(j)-1]
  readf, lun, line
  k = long(strsplit(line,delimiter,/extract))
  k = k[1:n_elements(k)-1]
  free_lun, lun

  if n_elements(j) eq 7 then begin
     readcol, filename, g,q0,q1,q2,q3,q4,q5,q6, format='d,d,d,d,d,d,d', skipline=2, delimiter=delimiter,/silent
     q_jk = [[q0],[q1],[q2],[q3],[q4],[q5],[q6]]
  endif else if n_elements(j) eq 8 then begin
     readcol, filename, g,q0,q1,q2,q3,q4,q5,q6,q7, format='d,d,d,d,d,d,d,d', skipline=2, delimiter=delimiter,/silent
     q_jk = [[q0],[q1],[q2],[q3],[q4],[q5],[q6],[q7]]
  endif else begin
     message, 'ERROR in expected coefficient file format'
  endelse
  sz = size(q_jk)
  n=sz[1]
  m=sz[2]

end

;; get the zero point correction to the Gaia EDR3 parallax 
;; (with minimal error checking)
function get_zpt, phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolor, ecl_lat, astrometric_params_solved

  path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir='sed') + path_sep()
  if astrometric_params_solved eq 31 then begin
     ;; 5 parameter solution
     color = nu_eff_used_in_astrometry
     read_coeffs,path+'z5_200720.txt', delimiter=' ', j=j, k=k, g=g, q_jk=q_jk, n=n, m=m
  endif else if astrometric_params_solved eq 95 then begin
     ;; 6 parameter solution
     color = pseudocolor 
     read_coeffs,path+'z6_200720.txt', delimiter=',', j=j, k=k, g=g, q_jk=q_jk, n=n, m=m
  endif

  sinbeta = sin(ecl_lat*!dpi/180d0)

  ;; basis functions evaluated at color and ecl_lat
  c = [1d0,$
       max([-0.24d0, min([0.24d0, color - 1.48d0])]),$
       min([0.24d0, max([0d0, 1.48d0 - color])])^3,$
       min([0d0, color - 1.24d0]),$
       max([0d0, color - 1.72d0])]

  b = [1d0, sinBeta, sinBeta^2 - 1d0/3d0]

  ;; coefficients must be interpolated between g(left) and g(left+1)
  ;; find the bin in g to the left of phot_g_mean_mag
  if phot_g_mean_mag le 6.0 then ig = 0 $
  else if phot_g_mean_mag gt 20 then ig = 11 $
  else begin
     allbigger = (where(phot_g_mean_mag ge g))
     ig = allbigger[n_elements(allbigger)-1]
  endelse

  ;; interpolate coefficients to gMag:
  h = max([0d0, min([1d0, (phot_g_mean_mag - g[ig]) / (g[ig + 1] - g[ig])])])

  ;; sum over the product of the coefficients to get the zero-point
  zpt = 0d0
  for i=0L, m-1 do zpt += ((1d0 - h)*q_jk[ig, i] + h*(q_jk[ig + 1, i]))*c[j[i]]*b[k[i]]

  return, zpt/1d3

end

pro getpars, ticid, dist=dist

;; TOI-201
ticid = '350618622'

if n_elements(dist) eq 0 then dist = 60d0

qtic = Exofast_Queryvizier('IV/38/tic','TIC ' + strtrim(ticid,2),/allcolumns,/cfa)
if (size(qtic))[2] ne 8 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

match = (where(qtic.tic eq ticid))[0]
if match eq -1 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

qtic = qtic[match]
star = [qtic.raj2000,qtic.dej2000]
gaiaid = qtic.gaia

;; EDR3
qgaia3=Exofast_Queryvizier('I/350/gaiaedr3',star,dist/60.,/silent,cfa=cfa,/all)

if (size(qgaia3))[2] eq 8 then begin
   match = (where(qgaia3.source eq gaiaid))[0]
   if match ne -1 then begin
      qgaia3 = qgaia3[match]

      phot_g_mean_mag = qgaia3.gmag 
      nu_eff_used_in_astrometry = qgaia3.nueff
      pseudocolor = qgaia3.pscol
      ecl_lat = qgaia3.elat
      astrometric_params_solved = qgaia3.solved

print, phot_g_mean_mag
print, nu_eff_used_in_astrometry
print, pseudocolor
print, ecl_lat
print, astrometric_params_solved



   endif
endif


end

pro test_all

  minmag = 6d0
  maxmag = 21d0
  nsteps = 10

  phot_g_mean_mag = minmag + (maxmag-minmag)*(dindgen(nsteps)/(nsteps-1))
  astrometric_params_solved = [31,95]

;  phot_g_mean_mag = 8.9492890d0
  nu_eff_used_in_astrometry = 1.59d0
  pseudocolor = 1.5d0
  ecl_lat = -78.283186d0
;  astrometric_params_solved = 31
;  astrometric_params_solved = 95
  


  for i=0L, n_elements(astrometric_params_solved)-1 do begin
     for j=0L, n_elements(phot_g_mean_mag)-1 do begin
        zpt = get_zpt(phot_g_mean_mag[j], nu_eff_used_in_astrometry, pseudocolor, ecl_lat, astrometric_params_solved[i])
        print, "gmag = " + strtrim(phot_g_mean_mag[j],2) + ' soltype=' + strtrim(astrometric_params_solved[i],2) + ': ' + strtrim(zpt,2)
     endfor
  endfor

end

