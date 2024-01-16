pro fitkepler, ticid
ticid = '399860444'

readcol, 'tickic_lookup.dat', ticid, kicid, format='a,a'
shuffle = sort(randomu(seed,n_elements(ticid)))
priormap = [['EEP','eep'],$
            ['Teff','teff'],$
            ['Mass','mstar'],$
            ['parallax','parallax'],$
            ['[Fe/H]','feh'],$
            ['tc','tc'],$
            ['period','period'],$
            ['i','ideg'],$
            ['p','p']]

for i=0L, n_elements(ticid)-1 do begin
   thisticid = ticid[shuffle[i]]
   thiskicid = ticid[shuffle[i]]

   ;; generate the prior file from Phill's file
   readcol, 'TIC' + ticid + '_sviMS_vX.priors', parname, value, err, format='a,d,d', comment='#'
   priorfilename = 'TIC' + ticid + '.priors'

   openw, priorlun, priorfilename
   nplanets = 0L
   for j=0L, n_elements(parname)-1 do begin
      match = where((strsplit(parname[0,j],'_',/extract))[0] eq priormap[0,*])
      if match[0] ne -1 then begin
         
         if strpos(parname[0,j],'_') ne -1 then begin
            planetnum = (strsplit(parname[0,j],'_',/extract))[1]  
            ndx = '_'+planetnum
            if nplanets lt long(planetnum+1) then nplanets = long(planetnum+1)
         endif else ndx = ''
         if priormap[1,match[0]] eq 'parallax' then begin
            printf, priorlun, priormap[1,match[0]]+ndx, value[j], err[j]            
         endif else if priormap[1,match[0]] eq 'rstar' then begin
            printf, priorlun, priormap[1,match[0]]+ndx, 10^value[j]
         endif else begin
            printf, priorlun, priormap[1,match[0]]+ndx, value[j]
         endelse
      endif
   endfor
   free_lun, lun

   ;; Phill's sed file
   sedfilename = 'TIC' + ticid + '.sed'

   ;; grab the Kepler data
   cmd = 'scp minerva@minerva-nas:/volume1/data/kepler/kplr*' + kicid + '*llc.fits .'
   spawn, cmd
   fitsfiles = file_search('*' + kicid + 'llc.fits',count=nfiles)
   for j=0L, nfiles-1 do begin
      data = mrdfits(fitsfiles[i],1)

      time = data.time + 2454833d0 ;; BJD_TDB
      flux = data.pdcsap_flux
      fluxerr = data.pdcsap_flux_err
      
      ;; remove flagged bad data points
      good = where(finite(time) and time ne 0 and finite(flux) and data.sap_quality eq 0L)
      time = time[good]
      flux = flux[good]
      fluxerr = fluxerr[good]
   
      ;; trim just the transits
      phase = dblarr(nplanets,n_elements(time))
      use = []
      for k=0L, nplanets-1 do begin
         phase[k,*] = (time-tc[k]) mod period[k]
         toohigh = where(phase[k,*] gt period[k]/2d0)
         if toohigh[0] ne -1 then phase[k,toohigh] -= period[k]
         toolow = where(phase[k,*] le -period[k]/2d0)
         if toolow[0] ne -1 then phase[k,toolow] += period[k]
         use = [use,where(phase[k,*] gt (-2d0*duration[k]) and phase[k,*] lt (2d0*duration[k]))]            
      endfor
      ;; remove overlapping elements
      sorted = sort(use)
      use = use[sorted]
      use = use[uniq(use)]
      
      ;; flatten light curve
            norm = keplerspline(time, flux, breakp=breakp, ndays=0.75)
            forprint, (time)[use],(flux/norm)[use],(fluxerr/norm)[use],format='(f0.6,x,f0.6,x,f0.6)',textout=tranfile,/nocomment,/silent
        




   endfor

endif




kicid = '001429092'
fitsfiles = file_search('*' + kicid + 'llc.fits',count=nfiles)
for i=0L, nfiles-1 do begin
   data = mrdfits(fitsfiles[i],1)
endfor

summaryfile = 'toi-2018-10-12.csv'
summaryfile = 'toi-plus-2018-11-03.csv'
summaryfile = 'toi-plus-2018-11-15.csv'
summaryfile = 'toi-plus-2018-11-30.csv'
summaryfile = 'toi-plus-2018-12-12.csv'
summaryfile = 'toi-plus-2018-12-27.csv'
summaryfile = 'toi-plus-2019-03-28.csv'
summaryfile = 'toi-plus-2019-04-24.csv'
maxsectors = 8

;tic = read_csv(summaryfile,header=tagnames, n_table_header=0)
;ticidndx = (where(tagnames eq 'tic_id'))[0]
;toiidndx = (where(tagnames eq 'toi_id'))[0]
;teffndx = (where(tagnames eq 'T_eff'))[0]
;uteffndx = (where(tagnames eq 'T_eff Err'))[0]
;srcndx = (where(tagnames eq 'src'))[0]
;depthndx = (where(tagnames eq 'Transit Depth'))[0]
;randx = (where(tagnames eq 'RA'))[0]
;decndx = (where(tagnames eq 'Dec'))[0]
;periodndx = (where(tagnames eq 'Period'))[0]
;nplanetsndx = (where(tagnames eq 'N_planets'))[0]
;rstarndx = (where(tagnames eq 'R_s'))[0]
;tcndx = (where(tagnames eq 'Epoc'))[0] ;; BJD_TDB - ???
;tmagndx = (where(tagnames eq 'Tmag'))[0]
;durationndx = (where(tagnames eq 'Duration'))[0]

;; this dictates the quality of the fit and runtime
maxsteps=7500
nthin=30
maxtime=2.5d0*86400d0

;; date string for the tgz files
caldat, systime(/julian,/utc), month,day,year
datestr = string(year,month, day,format='(i4,i02,i02)') 

;; get the files that get the lightcurves
for i=1L, maxsectors do begin
   filename = 'tesscurl_sector_' + strtrim(i,2) + '_lc.sh'
   if not file_test(filename) then begin
      spawn, 'wget http://archive.stsci.edu/missions/tess/download_scripts/sector/' + filename
   endif
endfor
   
spawn, 'hostname', output
hostname = (strsplit(output[0],'.',/extract))[0]

if hostname eq 'onion' then dir = '.' + path_sep() $
else dir = '/n/jjohnson_lab/jeastman/tess/' ;; slow, 4 TB, permanent
;   dir = '/n/regal/jjohnson_lab/' ;; scratch space -- fast, big, but deleted after 90 days

;; shuffle the files for super computer to prevent clashing threads
shuffle = sort(randomu(seed, n_elements(tic.(0))))

for i=0, n_elements(tic.(0))-1 do begin

   if n_elements(dothisticid) ne 0 then i = (where(shuffle eq (where(tic.(ticidndx) eq dothisticid))[0]))[0]

   ;; if it's not the first planet in the system, skip it 
   ;; we'll model it globally with the others 
   toiid = strtrim(tic.(toiidndx)[shuffle[i]],2)
   if strpos(toiid, '.01') eq -1 then continue

   ticid = strtrim(tic.(ticidndx)[shuffle[i]],2)
   path = dir + ticid + path_sep()
   prefix = path + ticid + '.'

;;  if it's from the quick look pipeline, skip it
;   if tic.(srcndx)[shuffle[i]] ne 'qlp' then stop ;continue

   ;; find the corresponding light curve files (postage stamps only?)
   spawn, 'grep ' + ticid + ' tesscurl_sector_*_lc.sh', output
   for j=0L, n_elements(output)-1 do begin
      entries = strsplit(output[j],'/ ',/extract)
      filename = entries(n_elements(entries)-1)
      if not file_test(dir + filename) then begin
         url = (strsplit(output[j],/extract))[6]
         spawn, 'wget -O ' + dir + filename + ' "' + url + '"'
      endif
   endfor
   
   fitsfiles = file_search(dir + 'kplr' + kicid + '*_llc.fits', count=nfiles)
   if nfiles eq 0 then begin
      print, 'No light curves for ' + kicid + ' (' + tic.(srcndx)[shuffle[i]] + ')'
      continue
   endif

   ;; if another thread is doing this one, skip it
   spawn, 'squeue -u jeastman', output
   for j=1L, n_elements(output)-1 do begin
      jobhostname = (strsplit(output[j],/extract))[7]
      if file_test(ticid + '_' + jobhostname + '.running') then goto, next
   endfor
   
   ;; if it's not running on any other threads, 
   ;; clean up files from pre-empted runs
   ;oldfiles = file_search(ticid + '_*.running', count=nold)
   ;if nold gt 0 then file_delete, oldfiles, /allow_nonexistent

   ;; create a file to tell other threads this one is doing it
   ;; do it here instead of after to reduce race conditions
   runningfile = ticid + '_' + hostname + '.running'
   spawn, 'touch ' + runningfile

   ;; if another thread has done this one, skip it
   if file_test(path) then begin

      ;; if it has run to completion and there is no additional data, skip it
      texfile = path + ticid + '.median.tex'
      logfile = path + ticid + '.log'

      ;; if the tex file exists, skip it
      if file_test(texfile) then begin
         ;; unless we have additional sectors, then add them and redo it
         junk = file_search(dir + '*.dat', count=nlcs)
         if nlcs eq nfiles and n_elements(dothisticid) eq 0 then begin
            file_delete, runningfile ;; delete this so it's not running
            continue
         endif else begin
            file_delete, texfile
            file_delete, logfile
         endelse
      endif      

      ;; if the log file has been updated in the last 4 hours, another
      ;; thread is probably still working on it; skip it
      ;if systime(/seconds)-file_modtime(logfile) lt 12d0*3600d0 and n_elements(dothisticid) eq 0 then continue
      
      ;; otherwise, redo the fit (previous attempt preempted?)
   endif else file_mkdir, path

   logfile = path + ticid + '.log'
   spawn, 'touch ' + logfile

   longcadence = tic.(srcndx)[shuffle[i]] eq 'qlp'

   ra = tic.(randx)[shuffle[i]]
   dec = tic.(decndx)[shuffle[i]]

   ;; find the other candidates in the system
   match = where(tic.(ticidndx) eq ticid, nplanets)
   p = sqrt(tic.(depthndx)[match]/1d6)
   tc = tic.(tcndx)[match] + 2457000d0
   period = tic.(periodndx)[match]
   duration = tic.(durationndx)[match]/24d0
   constants = mkconstants()
   ticids = (tic.(toiidndx))[match]
   
   priorfile = prefix + 'priors'
   tranfilebase = 'n20180701.TESS.TESS.' + ticid 
   sedfile = prefix + 'sed'

   ;if not file_test(tranfile + '.orig') then file_move, tranfile, tranfile+'.orig'
   for j=0L, nfiles-1 do begin

      fileparts = (strsplit(fitsfiles[j],'-_',/extract))
      if strpos(fileparts[2],'s') eq 0 and strlen(fileparts[2]) eq 5 then begin
         sector = 's' + string(long(strmid(fileparts[2],1,strlen(fileparts[2])-1)),format='(i02)')
      endif else if n_elements(fileparts) lt 9 then sector = 's01' $
      else sector = fileparts[8]

      tranfile =  path + tranfilebase + '.' + sector + '.dat'
      
      if not file_test(tranfile) then begin
      
         if tic.(srcndx)[shuffle[i]] eq 'qlp' then begin

            data = mrdfits(fitsfiles[j],1)

            ;; ***** generate the flattened light curve ****
            time = data.time + 2454833d0 
            flux = data.sap_flux
            
            ;; remove bad data points
            good = where(finite(time) and time ne 0 and finite(flux) and data.quality eq 0L)
            time = time[good]
            flux = flux[good]
            
            phase = dblarr(nplanets,n_elements(time))
            use = []
            for j=0L, nplanets-1 do begin
               phase[j,*] = (time-tc[j]) mod period[j]
               toohigh = where(phase[j,*] gt period[j]/2d0)
               if toohigh[0] ne -1 then phase[j,toohigh] -= period[j]
               toolow = where(phase[j,*] le -period[j]/2d0)
               if toolow[0] ne -1 then phase[j,toolow] += period[j]
               use = [use,where(phase[j,*] gt (-2d0*duration[j]) and phase[j,*] lt (2d0*duration[j]))]
            endfor
            
            ;; remove overlapping elements
            sorted = sort(use)
            use = use[sorted]
            use = use[uniq(use)]
            
            ;; everything not used can be used to calculate the error
            bad = time
            bad[use] = !values.d_nan
            oot = where(finite(bad))
            
            fluxerr = dblarr(n_elements(time)) + stdev(flux[oot])
            forprint, time[use],flux[use],fluxerr[use],format='(f0.6,x,f0.6,x,f0.6)',textout=tranfile,/nocomment,/silent
         endif else if tic.(srcndx)[shuffle[i]] eq 'spoc' then begin

            data = mrdfits(fitsfiles[j],1)
            time = data.time + 2454833d0 ;; BJD_TDB?
            flux = data.pdcsap_flux
            fluxerr = data.pdcsap_flux_err
            
            ;; remove flagged bad data points
            good = where(finite(time) and time ne 0 and finite(flux) and data.sap_quality eq 0L)
            time = time[good]
            flux = flux[good]
            fluxerr = fluxerr[good]
            
            ;; if the period is 0 (unknown), use all points
            unknownperiod = where(period eq 0d0)
            if unknownperiod[0] ne -1 then begin
               use = lindgen(n_elements(time)) 
               period[unknownperiod] = 30d0
            endif else begin
               phase = dblarr(nplanets,n_elements(time))
               use = []
               for k=0L, nplanets-1 do begin
                  phase[k,*] = (time-tc[k]) mod period[k]
                  toohigh = where(phase[k,*] gt period[k]/2d0)
                  if toohigh[0] ne -1 then phase[k,toohigh] -= period[k]
                  toolow = where(phase[k,*] le -period[k]/2d0)
                  if toolow[0] ne -1 then phase[k,toolow] += period[k]
                  use = [use,where(phase[k,*] gt (-2d0*duration[k]) and phase[k,*] lt (2d0*duration[k]))]            
               endfor
               ;; remove overlapping elements
               sorted = sort(use)
               use = use[sorted]
               use = use[uniq(use)]
            endelse               

            ;; flatten light curve
            norm = keplerspline(time, flux, breakp=breakp, ndays=0.75)
            forprint, (time)[use],(flux/norm)[use],(fluxerr/norm)[use],format='(f0.6,x,f0.6,x,f0.6)',textout=tranfile,/nocomment,/silent
            
         endif else begin
            print, 'source (' + strtrim(tic.(srcndx)[shuffle[i]],2) + ') not recognized; skipping'
            file_delete, runningfile
            continue
         endelse
      endif
   endfor

   ;; make the SED file (and grab the parallax while we're at it)
   if not file_test(sedfile) or not file_test(priorfile) then begin
      sedname = prefix + 'sed'
      mksed, tic.(ticidndx)[shuffle[i]], sedfile, ra=ra, dec=dec, parallax=parallax, uparallax=uparallax, dist=30,/over
   endif
   
   ;; ****** generate the prior file ********
   if not file_test(priorfile) then begin

      ;; query other parameters from the TIC
      mastobj = "{'service':'Mast.Catalogs.Tic.Cone','params':{'ra':" + strtrim(ra,2) + ",'dec':" + strtrim(dec,2) + ",'radius':0.001},'format':'json','timeout':10}"

      characters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
      randomstr = strjoin(characters[randomu(seed,12)*n_elements(characters)])

      filename = randomstr + '.txt'
      spawn, 'wget -O ' + filename + ' "https://mast.stsci.edu/api/v0/invoke?request=' + mastobj + '"'
      line = ''
      openr, lun, filename,/get_lun
      readf, lun, line
      free_lun, lun
      result = json_parse(line)
      rstar = ((result['data'])[0])['rad']
      urstar = ((result['data'])[0])['e_rad']
      mstar = ((result['data'])[0])['mass']
      umstar = ((result['data'])[0])['e_mass']

      contrastratio = ((result['data'])[0])['contratio']
      if n_elements(contrastratio) eq 1 then dilute = contrastratio/(1d0+contrastratio) $
      else dilute = 0d0
      feh = ((result['data'])[0])['MH']
      ufeh = ((result['data'])[0])['e_MH']
      teff = ((result['data'])[0])['Teff']
      uteff = ((result['data'])[0])['e_Teff']
      
      file_delete, filename

      eep = 355.65d0
      age = 4.6d0
      if n_elements(feh) eq 0 then mistfeh = 0d0 else mistfeh = feh
      if n_elements(mstar) eq 0 then mistmstar = 1d0 else mistmstar = mstar
      if n_elements(rstar) eq 0 then mistrstar = 1d0 else mistrstar = rstar
      if n_elements(teff) eq 0 then mistteff = 6000d0 else mistteff = teff
      mistchi2 = massradius_mist(eep,mistmstar,mistfeh,age,mistteff,mistrstar,mistfeh)
      ;; if we can't start here, start at solar (maybe just change eep?)
      if ~finite(mistchi2) then begin
         mstar = 1d0
         rstar = 1d0
      endif

      openw, lun, priorfile, /get_lun
      
      if n_elements(teff) ne 0 then begin
         ;if n_elements(uteff) ne 0 then printf, lun, 'teff', teff, uteff, format='(a,x,f0.4, x, f0.4)' $
         printf, lun, 'teff', teff, format='(a,x,f0.4)'
      endif
      if n_elements(feh) ne 0 then begin
         ;if n_elements(ufeh) ne 0 then printf, lun, 'feh', feh, ufeh, format='(a,x,f0.4, x, f0.4)' $
         printf, lun, 'feh', feh, format='(a,x,f0.4)'
      endif
;      if n_elements(rstar) ne 0 then printf, lun, 'rstar', rstar, format='(a,x,f0.4)'
;      if n_elements(mstar) ne 0 then printf, lun, 'mstar', mstar, format='(a,x,f0.4)'   
      if n_elements(dilute) ne 0 then printf, lun, 'dilute 0 ', dilute*0.1, format='(a,x,f0.8)' ;; dilution already accounted for in light curves, but allow 10% uncertainty
      
      junk = getavprior(ra=ra,dec=dec,line=line)
      printf, lun, line ;; AV prior
      
      if parallax gt 0d0 and finite(uparallax) then $
         printf, lun, 'parallax', parallax + 0.082d0, sqrt(uparallax^2 + 0.033d0^2), format='(a,x,f0.4, x, f0.4)'

      ;; the planet parameters from the DV reports
      if nplanets eq 0 then begin
         print, 'No planets found; skipping'
         file_delete, runningfile
         continue
      endif

      ticids = (tic.(toiidndx))[match]
      nplanets = 0
      for j=0, n_elements(ticids)-1 do nplanets = max([nplanets,double((strsplit(strtrim(ticids[j],2),'.',/extract))[1])/1000d0])

      ;nplanets = strtrim(tic.(nplanetsndx)[shuffle[i]],2)
      for j=0, nplanets-1 do begin
         printf, lun, string('tc_',j,tc[j],-1,tc[j]-period[j]/3d0,tc[j]+period[j]/3d0,format='(a,i1,x,f0.6,x,i2,x,f0.6,x,f0.6)')

         if period[j] ne 30d0 then printf, lun, string('period_',j,period[j],-1,period[j]*0.9,period[j]*1.1,format='(a,i1,x,f0.6,x,i2,x,f0.6,x,f0.6)') $
         else printf, lun, string('period_',j,period[j],format='(a,i1,x,f0.6)')
         printf, lun, 'p_' + strtrim(j,2), p[j], format='(a,x,f0.4)'
      endfor
      free_lun, lun
   ;; *********************************
   endif

   if longcadence then begin
      ninterp=dblarr(nfiles) + 10d0
      exptime=dblarr(nfiles) + 30d0

      if nfiles gt 1 then begin
         lctxt = ', ninterp=[' + string(lonarr(nfiles)+10,format='(' + strtrim(nfiles-1,2) + '(i2,",")' + ',i2,"]")') +$
                 ', exptime=[' + string(lonarr(nfiles)+30,format='(' + strtrim(nfiles-1,2) + '(i2,"d0,")' + ',i2,"d0]")')
      endif else lctxt = ', ninterp=[10],exptime=[30d0]'

   endif else begin
      ninterp=1d0
      exptime=2d0
      if nfiles gt 1 then begin
         lctxt = ', ninterp=[' + string(lonarr(nfiles)+1,format='(' + strtrim(nfiles-1,2) + '(i1,",")' + ',i1,"]")') +$
                 ', exptime=[' + string(lonarr(nfiles)+2,format='(' + strtrim(nfiles-1,2) + '(i1,"d0,")' + ',i1,"d0]")')
      endif else lctxt = ', ninterp=[1],exptime=[2d0]'
   endelse

   ;; make a pro file for future fits
   profilename = path + 'fit' + ticid + '.pro'
   if ~file_test(profilename) then begin
      openw, lun, profilename,/get_lun
      printf, lun, "pro fit" + ticid + ", debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin"
      printf, lun, ""
      printf, lun, "if n_elements(maxsteps) eq 0 then maxsteps=" + strtrim(maxsteps,2)
      printf, lun, "if n_elements(nthin) eq 0 then nthin=" + strtrim(nthin,2)
      printf, lun, ""
      printf, lun, "exofastv2, nplanets=" + strtrim(string(nplanets,format='(i)'),2) + ", tranpath='*.dat', priorfile='" + ticid + ".priors2',$"
      printf, lun, "           prefix='fitresults/" + ticid + ".', fluxfile='" + ticid + ".sed'" + lctxt + ", fitdilute=['TESS'],$"
      printf, lun, "           debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin"
      printf, lun, ""
      printf, lun, "end"
      free_lun, lun
   endif

   ;; ****  run the fit ****
   exofastv2, nplanets=nplanets, tranpath=path + tranfilebase + '.*.dat', priorfile = priorfile, $
              prefix=prefix, fluxfile=sedfile, ninterp=ninterp, exptime=exptime, fitdilute=['TESS'], $
              debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, $
              ntemp=8, maxtime=maxtime, rejectflatmodel=[lonarr(nfiles)+1]
   ;; *********************************

   ;; refine priors for next run
   mkprior, filename=prefix + 'mcmc.idl', priorfilename=priorfile + '2'
   
   ;; wrap up the results for upload to ExoFOP-TESS
   cd, path
;   spawn, 'pdflatex ' + ticid + '.median.tex' ;; compile latex table
   spawn, 'tar -cvzf TIC' + ticid + 'O-je' + datestr + '.tgz ' +$
          ticid + '.priors ' + $       ;; input priors
          ticid + '.priors2 ' + $      ;; refined priors
          'fit' + ticid + '.pro ' + $  ;; IDL file to refit planet
          ticid + '.sed ' +$           ;; input SED
          tranfilebase + '.*.dat ' + $ ;; input lightcurves
          ticid + '.mcmc.*ps ' + $     ;; final model plots
;          ticid + '.median.pdf ' + $   ;; pdf of summary table
          ticid + '.median.tex ' + $   ;; latex source for summary table
          ticid + '.median.csv ' + $   ;; machine readable summary table
          ticid + '.pdf.ps '     + $   ;; PDFs of all parameters
          ticid + '.covar.ps '   + $   ;; Corner plot of all fitted parameters
          ticid + '.chain.ps '   + $   ;; chains for all parameters
          ticid + '.mcmc.idl '   + $   ;; idl save file for all parameters
          ticid + '.log'               ;; log file for EXOFASTv2 run
   cd, '../'

   file_delete, runningfile ;; delete this so all threads know it's not running

   next:

endfor


end
