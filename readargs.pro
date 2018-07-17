;; Reads a file of input arguments to EXOFASTv2
;; Required for virtual machine (free) use
pro readargs, argfile, priorfile=priorfile, $
              rvpath=rvpath, tranpath=tranpath, dtpath=dtpath, fluxfile=fluxfile,$
              prefix=prefix,$
              circular=circular,fitslope=fitslope, secondary=secondary, $
              rossiter=rossiter,chen=chen,$
              fitthermal=fitthermal, fitreflect=fitreflect, fitdilute=fitdilute,$
              nthin=nthin, maxsteps=maxsteps, dontstop=dontstop, $
              debug=debug, stardebug=stardebug, verbose=verbose, randomfunc=randomfunc, seed=seed,$
              bestonly=bestonly, plotonly=plotonly, refinestar=refinestar, $
              longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
              maxgr=maxgr, mintz=mintz, $
              yy=yy, torres=torres, nomist=nomist, noclaret=noclaret, tides=tides, nplanets=nplanets, $
              fitrv=fitrv, fittran=fittran,fitdt=fitdt,lineark=lineark,$
              ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
              earth=earth, i180=i180, nocovar=nocovar,alloworbitcrossing=alloworbitcrossing,stretch=stretch

line = ''
openr, lun, argfile, /get_lun
while not eof(lun) do begin
   readf, lun, line
   line = (strsplit(line,'#',/extract,/preserve_null))[0]
   if line ne '' then begin
      entries = strsplit(line,'=',/extract)
      if n_elements(entries) eq 2 then begin

         ;; replace '$EXOFAST_PATH' with its value
         entries[1] = strjoin(strsplit(entries[1],'$EXOFAST_PATH',/regex,/extract,/preserve_null),getenv('EXOFAST_PATH'))

         if strupcase(strtrim(entries[0],2)) eq 'PRIORFILE' then begin
            priorfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'RVPATH' then begin
            rvpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'TRANPATH' then begin
            tranpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'DTPATH' then begin
            dtpath = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'FLUXFILE' then begin
            fluxfile = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'PREFIX' then begin
            prefix = strtrim(entries[1],2)
         endif else if strupcase(strtrim(entries[0],2)) eq 'CIRCULAR' then begin
            circular = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITSLOPE' then begin
            fitslope = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'SECONDARY' then begin
            secondary = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'ROSSITER' then begin
            rossiter = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'CHEN' then begin
            chen = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITTHERMAL' then begin
            fitthermal = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITREFLECT' then begin
            fitreflect = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITDILUTE' then begin
            fitdilute = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NTHIN' then begin
            nthin = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'MAXSTEPS' then begin
            maxsteps = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'DONTSTOP' then begin
            dontstop = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'DEBUG' then begin
            debug = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'STARDEBUG' then begin
            stardebug = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'VERBOSE' then begin
            verbose = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'SEED' then begin
            seed = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'BESTONLY' then begin
            bestonly = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'PLOTONLY' then begin
            plotonly = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'REFINESTAR' then begin
            refinestar = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'LONGCADENCE' then begin
            longcadence = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'EXPTIME' then begin
            exptime = double(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NINTERP' then begin
            ninterp = long(entries[1])
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
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOCLARET' then begin
            noclaret = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TIDES' then begin
            tides = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'NPLANETS' then begin
            nplanets = long(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITRV' then begin
            fitrv = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITTRAN' then begin
            fittran = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'FITDT' then begin
            fitdt = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'LINEARK' then begin
            lineark = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'TTVS' then begin
            ttvs = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TIVS' then begin
            tivs = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'TDELTAVS' then begin
            tdeltavs = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'EARTH' then begin
            earth = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'I180' then begin
            i180 = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'NOCOVAR' then begin
            nocovar = boolean(entries[1])
         endif else if strupcase(strtrim(entries[0],2)) eq 'ALLOWORBITCROSSING' then begin
            alloworbitcrossing = boolean(json_parse(entries[1],/toarray))
         endif else if strupcase(strtrim(entries[0],2)) eq 'STRETCH' then begin
            stretch = boolean(entries[1])
         endif else begin
            print, entries[0] + ' argument not recognized'
         endelse
      endif else print, 'invalid line: ' + line
   endif
endwhile
free_lun, lun

end
