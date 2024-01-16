;; Reads a file of input arguments to EXOFASTv2
;; Required for virtual machine (free) use
pro readargs, argfile, priorfile=priorfile, $
              prefix=prefix,$
              rvpath=rvpath, tranpath=tranpath, $
              astrompath=astrompath, dtpath=dtpath, $
              fluxfile=fluxfile,mistsedfile=mistsedfile,$
              sedfile=sedfile,specphotpath=specphotpath,$
              fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor, $
              fehsedfloor=fehsedfloor, oned=oned,$
              yy=yy, nomist=nomist, parsec=parsec, $
              torres=torres, mannrad=mannrad, mannmass=mannmass, $
              teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
              rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
              fitthermal=fitthermal, fitellip=fitellip, $
              fitreflect=fitreflect, fitphase=fitphase, $
              fitbeam=fitbeam, derivebeam=derivebeam, $
              nstars=nstars,starndx=starndx, $
              diluted=diluted,fitdilute=fitdilute, $
              nplanets=nplanets, $
              fittran=fittran, fitrv=fitrv, $
              rossiter=rossiter, fitdt=fitdt, $
              circular=circular, tides=tides, $ 
              alloworbitcrossing=alloworbitcrossing, $
              chen=chen, i180=i180, $
              fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
              noclaret=noclaret, $
              ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
              longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
              rejectflatmodel=rejectflatmodel,$
              fitspline=fitspline, splinespace=splinespace, $
              fitwavelet=fitwavelet, $              
              fitlogmp=fitlogmp,$
              novcve=novcve, nochord=nochord, fitsign=fitsign, $
              fittt=fittt, earth=earth, $             
              transitrange=transitrange,rvrange=rvrange,$
              sedrange=sedrange,emrange=emrange, $
              debug=debug, verbose=verbose, delay=delay,$
              maxsteps=maxsteps, nthin=nthin, maxtime=maxtime, $
              maxgr=maxgr, mintz=mintz, $
              dontstop=dontstop, $              
              ntemps=ntemps,tf=tf,keephot=keephot,$
              randomfunc=randomfunc, seed=seed,$
              stretch=stretch,$
              skiptt=skiptt, $
              usernote=usernote,$
              mksummarypg=mksummarypg,$
              nocovar=nocovar,$
              plotonly=plotonly, bestonly=bestonly,$
              logname=logname

line = ''
openr, lun, argfile, /get_lun
while not eof(lun) do begin
   readf, lun, line
   line = (strsplit(line,'#',/extract,/preserve_null))[0]
   if line ne '' then begin
      entries = strsplit(line,'=',/extract)
      if n_elements(entries) eq 2 then begin

         ;; replace environment variables with their values
         entries[1] = exofast_expand_path(entries[1])

         if strupcase(strtrim(entries[0],2)) eq 'PRIORFILE' then begin
            priorfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'RVPATH' then begin
            rvpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'TRANPATH' then begin
            tranpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'ASTROMPATH' then begin
            astrompath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'DTPATH' then begin
            dtpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FLUXFILE' then begin
            fluxfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'MISTSEDFILE' then begin
            mistsedfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'SEDFILE' then begin
            sedfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'SPECPHOTPATH' then begin
            specphotpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FBOLSEDFLOOR' then begin
            fbolsedfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TEFFSEDFLOOR' then begin
            teffsedfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FEHSEDFLOOR' then begin
            fehsedfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'ONED' then begin
            oned = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TEFFEMFLOOR' then begin
            teffemfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FEHEMFLOOR' then begin
            fehemfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'RSTAREMFLOOR' then begin
            rstaremfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'AGEEMFLOOR' then begin
            ageemfloor = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'PREFIX' then begin
            prefix = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'CIRCULAR' then begin
            circular = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITSLOPE' then begin
            fitslope = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITQUAD' then begin
            fitquad = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITTT' then begin
            fittt = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'RVEPOCH' then begin
            rvepoch = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'ROSSITER' then begin
            rossiter = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'CHEN' then begin
            chen = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'DILUTED' then begin
            diluted = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITDILUTE' then begin
            fitdilute = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITTHERMAL' then begin
            fitthermal = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITREFLECT' then begin 
            fitreflect = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITPHASE' then begin 
            fitphase = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITELLIP' then begin
            fitellip = json_parse(entries[1],/toarray)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITBEAM' then begin
            fitbeam = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'DERIVEBEAM' then begin
            derivebeam = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NTHIN' then begin
            nthin = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MAXSTEPS' then begin
            maxsteps = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MAXTIME' then begin
            maxtime = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'DONTSTOP' then begin
            dontstop = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NTEMPS' then begin
            ntemps = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TF' then begin
            tf = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'KEEPHOT' then begin
            keephot = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'DEBUG' then begin
            debug = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'VERBOSE' then begin
            verbose = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'DELAY' then begin
            delay = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'SEED' then begin
            seed = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'PLOTONLY' then begin
            plotonly = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'BESTONLY' then begin
            bestonly = boolean(entries[1])            
         endif else if strupcase(strtrim(entries[0],2)) eq 'LONGCADENCE' then begin
            longcadence = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'EXPTIME' then begin
            exptime = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NINTERP' then begin
            ninterp = long(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'REJECTFLATMODEL' then begin
            rejectflatmodel = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MAXGR' then begin
            maxgr = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MINTZ' then begin
            mintz = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'YY' then begin
            yy = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TORRES' then begin
            torres = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOMIST' then begin
            nomist = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'PARSEC' then begin
            parsec = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MANNRAD' then begin
            mannrad = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MANNMASS' then begin
            mannmass = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOCLARET' then begin
            noclaret = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'TIDES' then begin
            tides = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NPLANETS' then begin
            nplanets = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NSTARS' then begin
            nstars = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'STARNDX' then begin
            starndx = long(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITRV' then begin
            fitrv = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITTRAN' then begin
            fittran = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITDT' then begin
            fitdt = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITLOGMP' then begin
            fitlogmp = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'TTVS' then begin
            ttvs = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'TIVS' then begin
            tivs = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'TDELTAVS' then begin
            tdeltavs = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'EARTH' then begin
            earth = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TRANSITRANGE' then begin
            transitrange = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'RVRANGE' then begin
            rvrange = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'SEDRANGE' then begin
            sedrange = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'EMRANGE' then begin
            emrange = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'I180' then begin
            i180 = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOVCVE' then begin
            novcve = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOCHORD' then begin
            nochord = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITSIGN' then begin
            fitsign = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOCOVAR' then begin
            nocovar = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'ALLOWORBITCROSSING' then begin
            alloworbitcrossing = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'STRETCH' then begin
            stretch = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITSPLINE' then begin
            fitspline = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'SPLINESPACE' then begin
            splinespace = double(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITWAVELET' then begin
            fitwavelet = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'SKIPTT' then begin
            skiptt = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'USERNOTE' then begin
            usernote = entries[1]
         endif else if strupcase(strtrim(entries[0],2)) eq 'MKSUMMARYPG' then begin
            mksummarypg = boolean(entries[1])
         endif else begin
            printandlog, entries[0] + ' argument not recognized', logname
         endelse
      endif else printandlog, 'invalid line: ' + line, logname
   endif
endwhile
free_lun, lun

end
