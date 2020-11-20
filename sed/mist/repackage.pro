pro repackage

;; unpack all the files 
;; downloaded from http://waps.cfa.harvard.edu/MIST/model_grids.html
;txzfiles = file_search('*.txz',count=nfiles)
;for i=0L, nfiles-1 do spawn, 'xz -d < ' + txzfiles[i] + ' | tar xvf -'

systems=['CFHT','DECam','GALEX','HST_ACSHR','HST_ACSWF','HST_WFC3',$
        'HST_WFPC2','JWST','LSST','PanSTARRS','SDSS','SPITZER',$
        'SkyMapper','Swift','UBVRIplus','UKIDSS','WFIRST','WISE',$
        'WashDDOuvby']
systems=['UBVRIplus']

;; meta data for each filter
readcol, 'Filters.csv', $
         fname2, lambda_mean, lambda_eff, lambda_min, lambda_max, w_eff, zp, $
         delimiter=',',skipline=1, format='a,d,d,d,d,d,d'
for i=0L, n_elements(fname2)-1 do fname2[i] = (strsplit(fname2[i],'/',/extract))[-1]

nsystems = n_elements(systems)

nteff = 70
nlogg = 26
nfeh = 18
nav = 13

for n=0L, nsystems-1 do begin
   
   hdr = ''
   bcfiles = file_search('feh*.' + systems[n] + '*',count=nfiles)

   ;; meta data containing MIST assumed zero points and mag types
   readcol, systems[n], junk, fname, type, vega_zp, st_zp, ab_zp, vegast, vegaab, format='a,a,a,d,d,d,d,d';, skipline=4
   
   for i=0L, nfiles-1 do begin


      ;; get the number of filters
      if i eq 0L then begin
         openr, lun, bcfiles[i], /get_lun
         for junk=0L, 5 do readf, lun, hdr
         hdr = strsplit(hdr,/extract)
         free_lun, lun
         nfilters = n_elements(hdr)-6
         array = dblarr(nteff, nlogg, nfeh, nav, nfilters)         
      endif

      ;; remove comments, read file
      spawn, 'grep -v "#" ' + bcfiles[i], output
      nlines = 23660L
      teff = dblarr(nlines)
      logg = dblarr(nlines)
      feh = dblarr(nlines)
      av = dblarr(nlines)
      rv = dblarr(nlines)
      filters = dblarr(nlines,nfilters)
      for j=0L, n_elements(output)-1 do begin
         entries = double(strsplit(output[j],/extract))
         teff[j] = entries[0]
         logg[j] = entries[1]
         feh[j] = entries[2]
         av[j] = entries[3]
         rv[j] = entries[4]
         filters[j,*] = entries[5:*]
      endfor

      teffsort = teff[sort(teff)]
      uniqteffs = teffsort[uniq(teffsort)]
      
      loggsort = logg[sort(logg)]
      uniqloggs = loggsort[uniq(loggsort)]
      
      avsort = av[sort(av)]
      uniqavs = avsort[uniq(avsort)]

;teffgrid = uniqteffs
;logggrid = uniqloggs
;avgrid = uniqavs
;fehgrid = [-4d0,-3.5d0,[-3d0 + dindgen(16)*0.25d0]]
;      save, teffgrid, logggrid, fehgrid, avgrid, filename='mist.sed.grid.idl'
;
;stop
      
      ;; populate the 4D array
      for j=0L, nteff-1 do begin
         for k=0L, nlogg-1 do begin
            for l=0L, nav-1 do begin
               match = where(teff eq uniqteffs[j] and logg eq uniqloggs[k] and av eq uniqavs[l],nmatch)
               if nmatch ne 1 then stop
               for m=0L, nfilters-1 do begin
                  array[j,k,i,l,m] = filters[match,m]
                  array[j,k,i,l,m] = filters[match,m]
               endfor
            endfor
         endfor
      endfor
      
   endfor

   for m=0L, nfilters-1 do begin
      bcarray = array[*,*,*,*,m]
      match = where(fname eq hdr[6+m])
      match2 = where(fname2 eq hdr[6+m])

      if match[0] eq -1 or match2[0] eq -1 then begin
         print, match, match2, hdr[6+m]
         stop
      endif

      filterproperties = create_struct('name',hdr[6+m],$
                                       'lambda_eff',lambda_eff[match2],$
                                       'w_eff',w_eff[match2],$
                                       'type',type[match],$
                                       'vega_zp',vega_zp[match],$
                                       'st_zp',st_zp[match],$
                                       'ab_zp',ab_zp[match],$
                                       'vegast',vegast[match],$
                                       'vegaab',vegaab[match])
      save, bcarray, filterproperties,filename=hdr[6+m] + '.idl'
      print, 'Done with ' + hdr[6+m] + '.idl'
   endfor

   print, 'Done with ' + systems[n]
   
endfor

end
