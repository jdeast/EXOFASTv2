;+
; NAME:
;   MKSS
;
; PURPOSE:
; Creates a stellar system structure that describes an arbitrary
; number of stars, planets, observed bands, rv telescopes (a new zero
; point for each), and observed transits.  This structure is the input
; to basically everything in EXOFASTv2 and is designed to be easily
; extensible. To add a new derived parameter to the output table, add
; it here, then calculate it in DERIVEPARS. To add a new fitted
; parameter, add it here, then use it in EXOFAST_CHI2V2 to influence
; the model or the likelihood/chi^2 directly.
;
; Priors may be applied to any parameter (fitted or derived).
;
; CALLING SEQUENCE:
;
;   This code is deep in the bowels of EXOFASTv2 and is not intended
;   to be user callable.
;
; INPUTS:
;  See $EXOFAST_PATH/exofastv2.pro for an explanation of most inputs.
; 
;  NVALUES  - By default, parameter.value is a scalar. Set this
;             to make it an array of length NVALUES.
;
;  SILENT   - If set, supress output messages.
;
;  CHI2FUNC - A string specifying the name of the function that
;             computes the model's chi^2.
;
; MODIFICATION HISTORY
; 
;  2023/12 -- Documentation cleanup.
;-
function mkss, priorfile=priorfile, $
               prefix=prefix,$
               ;; data file inputs
               rvpath=rvpath, tranpath=tranpath, $
               astrompath=astrompath, dtpath=dtpath, $
               ;; SED model inputs
               fluxfile=fluxfile, mistsedfile=mistsedfile, $
               sedfile=sedfile, specphotpath=specphotpath,$
               noavprior=noavprior,$
               fbolsedfloor=fbolsedfloor,teffsedfloor=teffsedfloor,$
               fehsedfloor=fehsedfloor, oned=oned,$
               ;; evolutionary model inputs
               yy=yy0, nomist=nomist, parsec=parsec0, $
               torres=torres0, mannrad=mannrad0, mannmass=mannmass0,$
               teffemfloor=teffemfloor, fehemfloor=fehemfloor, $
               rstaremfloor=rstaremfloor,ageemfloor=ageemfloor,$
               ;; BEER model inputs
               fitthermal=fitthermal, fitellip=fitellip, $
               fitreflect=fitreflect, fitphase=fitphase,$
               fitbeam=fitbeam, derivebeam=derivebeam, $
               ;; star inputs
               nstars=nstars, starndx=starndx, $
               diluted=diluted, fitdilute=fitdilute, $
               ;; planet inputs
               nplanets=nplanets, $
               fittran=fittran,fitrv=fitrv,$
               rossiter=rossiter, fitdt=fitdt,$ 
               circular=circular, tides=tides, $
               alloworbitcrossing=alloworbitcrossing,$
               chen=chen, i180=i180,$
               ;; RV inputs
               fitslope=fitslope, fitquad=fitquad, rvepoch=rvepoch, $
               ;; transit inputs
               noclaret=noclaret,$
               ttvs=ttvs, tivs=tivs, tdeltavs=tdeltavs, $
               longcadence=longcadence, exptime=exptime, ninterp=ninterp, $
               rejectflatmodel=rejectflatmodel,$
               fitspline=fitspline, splinespace=splinespace, $
               fitwavelet=fitwavelet, $            
               ;; reparameterization inputs
               fitlogmp=fitlogmp,$
               novcve=novcve, nochord=nochord, fitsign=fitsign, $
               fittt=fittt, earth=earth, $
               ;; plotting inputs
               transitrange=transitrange,rvrange=rvrange,$
               sedrange=sedrange,emrange=emrange, $
               ;; debugging inputs
               debug=debug, verbose=verbose, delay=delay, $
               ;; internal inputs
               nvalues=nvalues,$ 
               silent=silent, $
               chi2func=chi2func, $
               logname=logname, $
               best=best

if n_elements(transitrange) eq 0 then transitrange=dblarr(6)+!values.d_nan
if n_elements(rvrange) eq 0 then rvrange=dblarr(6) + !values.d_nan
if n_elements(sedrange) eq 0 then sedrange=dblarr(6) + !values.d_nan
if n_elements(emrange) eq 0 then emrange=dblarr(4) + !values.d_nan

if n_elements(nstars) eq 0 then nstars = 1L

if nstars eq 1 or keyword_set(noavprior) then avprior = 0B $
else avprior = 1B

if nstars lt 1 then begin
   printandlog, 'NSTARS must be greater than 0', logname
endif


if n_elements(nomist) eq 0 then begin
   ;; default uses MIST for all stars
   mist = bytarr(nstars)+1B
endif else if n_elements(nomist) eq 1 then begin
   if keyword_set(nomist) then begin
      ;; if set as a keyword, disable for all stars
      mist = bytarr(nstars)
   endif else begin
      mist = bytarr(nstars)+1B
   endelse
endif else if n_elements(nomist) eq nstars then begin
   ;; if set as an NSTARS array, use array
   mist = ~nomist
endif else begin
   printandlog, 'NOMIST must have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(parsec0) eq 0 then begin
   ;; don't use PARSEC by default
   parsec = bytarr(nstars)
endif else if n_elements(parsec0) eq 1 then begin
   if keyword_set(parsec0) then begin
      ;; if set as a keyword, disable for all stars
      parsec = bytarr(nstars) + 1B
   endif else begin
      parsec = bytarr(nstars)
   endelse
endif else if n_elements(parsec) eq nstars then begin
   ;; if set as an NSTARS array, use array
   parsec = parsec0
endif else begin
   printandlog, 'PARSEC must have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(torres0) eq 0 then begin
   ;; don't use TORRES by default
   torres = bytarr(nstars)
endif else if n_elements(torres0) eq 1 then begin
   if keyword_set(torres0) then begin
      ;; if set as a keyword, enable for all stars
      torres = bytarr(nstars) + 1B
   endif else begin
      torres = bytarr(nstars)      
   endelse
endif else if n_elements(torres) eq nstars then begin
   ;; if set as an NSTARS array, use array
   torres = torres0
endif else begin
   printandlog, 'TORRES mist have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(mannrad0) eq 0 then begin
   ;; don't use MANN by default
   mannrad = bytarr(nstars)
endif else if n_elements(mannrad0) eq 1 then begin
   mannrad = bytarr(nstars) + keyword_set(mannrad0)
endif else if n_elements(mannrad) eq nstars then begin
   ;; if set as an NSTARS array, use array
   mannrad = mannrad0
endif else begin
   printandlog, 'MANNRAD mist have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(mannmass0) eq 0 then begin
   ;; don't use MANN by default
   mannmass = bytarr(nstars)
endif else if n_elements(mannmass0) eq 1 then begin
   mannmass = bytarr(nstars) + keyword_set(mannmass0)
endif else if n_elements(mannmass) eq nstars then begin
   ;; if set as an NSTARS array, use array
   mannmass = mannmass0
endif else begin
   printandlog, 'MANNMASS mist have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(yy0) eq 0 then begin
   ;; don't use YY by default
   yy = bytarr(nstars)
endif else if n_elements(yy0) eq 1 then begin
   if keyword_set(yy0) then begin
      ;; if set as a keyword, disable for all stars
      yy = bytarr(nstars) + 1B
   endif else begin
      yy = bytarr(nstars)
   endelse
endif else if n_elements(yy) eq nstars then begin
   ;; if set as an NSTARS array, use array
   yy = yy0
endif else begin
   printandlog, 'YY mist have 0, 1, or NSTARS elements', lognname
   return, -1
endelse

if n_elements(fbolsedfloor) eq 0 then fbolsedfloor = 0.024d0
if n_elements(teffsedfloor) eq 0 then teffsedfloor = 0.020d0
if n_elements(fehsedfloor) eq 0 then fehsedfloor = 0d0

if n_elements(teffemfloor) eq 0 then teffemfloor = -1d0
if n_elements(fehemfloor) eq 0 then fehemfloor = -1d0
if n_elements(rstaremfloor) eq 0 then rstaremfloor = -1d0
if n_elements(ageemfloor) eq 0 then ageemfloor = -1d0

if n_elements(delay) eq 0 then delay = 0L $
else begin
   if delay ne 0 then $
      printandlog, 'You have specified a delay!! This is only designed to test threading and will needlessly slow down the fit. Are you sure?', logname
endelse

if max(mist + yy + parsec + torres + (mannrad or mannmass)) gt 1 then begin
   printandlog, 'You are **STRONGLY** advised to disable all but one evolutionary model (they are not independent), but type ".con" to proceed', logname
   return, -1
endif

if n_elements(nplanets) eq 0 then nplanets = 1
if n_elements(circular) ne nplanets and n_elements(circular) gt 1 then begin
   printandlog, "CIRCULAR must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fittran) ne nplanets and n_elements(fittran) gt 1 then begin
   printandlog, "FITTRAN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fitrv) ne nplanets and n_elements(fitrv) gt 1 then begin
   printandlog, "FITRV must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(chen) ne nplanets and n_elements(chen) gt 1 then begin
   printandlog, "CHEN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(i180) ne nplanets and n_elements(i180) gt 1 then begin
   printandlog, "I180 must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(rossiter) ne nplanets and n_elements(rossiter) gt 1 then begin
   printandlog, "ROSSITER must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fitdt) ne nplanets and n_elements(fitdt) gt 1 then begin
   printandlog, "FITDT must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fitlogmp) ne nplanets and n_elements(fitlogmp) gt 1 then begin
   printandlog, "FITLOGMP must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(novcve) ne nplanets and n_elements(novcve) gt 1 then begin
   printandlog, "NOVCVE must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(nochord) ne nplanets and n_elements(nochord) gt 1 then begin
   printandlog, "NOCHORD must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fitsign) ne nplanets and n_elements(fitsign) gt 1 then begin
   printandlog, "FITSIGN must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fitbeam) ne nplanets and n_elements(fitbeam) gt 1 then begin
   printandlog, "FITBEAM must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(derivebeam) ne nplanets and n_elements(derivebeam) gt 1 then begin
   printandlog, "DERIVEBEAM must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif
if n_elements(fittt) ne nplanets and n_elements(fittt) gt 1 then begin
   printandlog, "FITTT must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif

if n_elements(starndx) eq 0 then starndx = lonarr(nplanets>1)

if n_elements(starndx) ne nplanets and nplanets gt 0 then begin
   printandlog, "STARNDX must have NPLANETS (" + strtrim(nplanets,2) + ") elements",logname
   return, -1
endif

if not keyword_set(longcadence) then longcadence=0B
if n_elements(fitthermal) eq 0 then fitthermal = ['']
if n_elements(fitreflect) eq 0 then fitreflect = ['']
if n_elements(fitphase) eq 0 then fitphase = ['']
if n_elements(fitellip) eq 0 then fitellip = ['']
if n_elements(tranpath) eq 0 then tranpath = ''
if n_elements(rvpath) eq 0 then rvpath = ''
if n_elements(astrompath) eq 0 then astrompath = ''
if n_elements(dtpath) eq 0 then dtpath = ''
if n_elements(fluxfile) eq 0 then fluxfile = ''
if n_elements(mistsedfile) eq 0 then mistsedfile = ''
if n_elements(sedfile) eq 0 then sedfile = ''
if n_elements(specphotpath) eq 0 then specphotpath = ''

if tranpath ne '' then begin
   tranfiles=file_search(tranpath,count=ntran)
   if ntran eq 0 then begin
      printandlog, "No transit files files found matching " + strtrim(tranpath,2) + "; please check TRANPATH", logname
      return, -1
   endif
endif else begin
   ntran = 0
   tranpath = ''
endelse

if n_elements(noclaret) eq 1 then begin
   noclaret = bytarr(ntran) + keyword_set(noclaret)
endif else if n_elements(noclaret) ne ntran and n_elements(noclaret) ne 0 and ntran gt 0 then begin
   printandlog, 'NOCLARET has ' + strtrim(n_elements(noclaret),2) + ' elements; must be an NTRANSITS (' + strtrim(ntran,2) + ') element array', logname
   return, -1
end
if n_elements(noclaret) eq 0 then begin
   if ntran gt 0 then noclaret = bytarr(ntran) $
   else noclaret = [0B]
endif

if specphotpath ne '' then begin
   specfiles = file_search(specphotpath,count=nspecfiles)
   if nspecfiles eq 0 then begin          
      printandlog, "No spectrophotometry files files found matching " + strtrim(specphotpath,2) + "; please check SPECPHOTPATH", logname
      return, -1
   endif
   if not file_test(sedfile) then begin
      printandlog, 'Spectrophotometry is only supported with the NextGen SED model. Use SEDFILE to specify your SED file or remove SPECPHOTPATH'
      return, -1
   endif
endif else nspecfiles = 0

;; read in the transit files
if tranpath ne '' or astrompath ne '' then begin

   if astrompath ne '' then begin
      astromfiles=file_search(astrompath,count=nastrom)
      if nastrom eq 0 then begin
         printandlog, "No astrometry files files found matching " + strtrim(astrompath,2) + "; please check ASTROMPATH", logname
         return, -1
      endif
   endif else begin
      nastrom = 0
      astrompath = ''
   endelse

   ;; find the unique bands
   allowedbands = ['U','B','V','R','I','J','H','K',$
                   'Sloanu','Sloang','Sloanr','Sloani','Sloanz',$
                   'Kepler','TESS','CoRoT','Spit36','Spit45','Spit58','Spit80',$
                   'u','b','v','y']
   prettybands = ['U','B','V','R','I','J','H','K',$
                  "u'","g'","r'","i'","z'",$
                  'Kepler','TESS','CoRoT','$3.6\mu m$','$4.5\mu m$','$5.8\mu m$','$8.0\mu m$',$
                  'u','b','v','y']
   bands = strarr(ntran+nastrom)
   for i=0, ntran+nastrom-1 do begin
      if i lt ntran then begin
         bands[i] = (strsplit(file_basename(tranfiles[i]),'.',/extract))(1)
         checkband = ~keyword_set(noclaret[i])
      endif else begin
         bands[i] = (strsplit(file_basename(astromfiles[i-ntran]),'.',/extract))(1)
         checkband = 0
      endelse

      if bands[i] eq 'u' then begin
         bands[i] = 'Sloanu'
      endif else if bands[i] eq 'g' then begin
         bands[i] = 'Sloang'
      endif else if bands[i] eq 'r' then begin
         bands[i] = 'Sloanr'
      endif else if bands[i] eq 'i' then begin
         bands[i] = 'Sloani'
      endif else if bands[i] eq 'z' then begin
         bands[i] = 'Sloanz'
      endif
      if (where(allowedbands eq bands[i]))[0] eq -1 and checkband then begin
         printandlog, 'ERROR: band (' + bands[i] + ') not recognized from the transit file "' + tranfiles[i] + '"',logname
         printandlog, 'When CLARET tables are used for a transit, you must select from an allowed band',logname 
         printandlog, 'Please check the documentation for required filename format and select one of the following:',logname
         printandlog, string(allowedbands),logname
         return, -1
      endif
   endfor
   bands = bands[uniq(bands, sort(bands))]
   nband = n_elements(bands)
   if ~keyword_set(silent) then begin
      printandlog, 'The index for each fitted band is', logname
      for i=0L, nband-1 do printandlog, string(i,bands[i],format='(i2,x,a)'),logname
      printandlog, '', logname
   endif
endif else begin
   ntran = 0
   tranpath = ''
   nastrom = 0
   astrompath = ''
   nband = 0
endelse

;; initialize the diluted (NTRANxNSTARs) byte array
;; if not specified, default to no dilution
if n_elements(diluted) eq 0 then begin
   diluted = bytarr(ntran>1, nstars)
endif else if n_elements(diluted) eq 1 then begin
   ;; if specified as a keyword, dilute all LCs according to that keyword
   diluted = bytarr(ntran>1, nstars)+keyword_set(diluted)
endif

;; make sure the input make sense
sz = size(diluted) 
if total(diluted) gt 0 then begin
   if sz[0] eq 1 and nstars eq 1 and sz[2] eq ntran then begin  
      ;; only 1 star, cannot compute dilution
      printandlog, 'DILUTED cannot be used with only one star.',logname
      printandlog, 'Either model multiple stars to model dilution according to their SEDs', logname
      printandlog, 'or use FITDILUTE to FIT dilution independently for each lightcurve.', logname
      return, -1
   endif else if sz[0] eq 2 and ntran eq sz[1] and nstars eq sz[2] then begin
      ;; good (multiple stars), do nothing
   endif else begin
      printandlog, 'DILUTED (' + strtrim(sz[2],2) + ' x ' + strtrim(sz[1],2) + ') must be an NTRANxNSTARS (' + strtrim(ntran,2) + ' x ' + strtrim(nstars,2) + ') array', logname
      return, -1
   endelse
endif
dilutestarndx = where(total(diluted,1)) ;; this specifies which stars should be computed for deblending

if n_elements(fitdilute) eq 0 then begin
   ;; if not specified, we want to fit the dilution for every lightcurve specified in the DILUTED array
   ;; which is none, if DILUTED is not set)
   if nstars eq 1 then fitdilute = bytarr(ntran > 1) $
   else fitdilute = total(diluted,2)
endif else if n_elements(fitdilute) eq 1 then begin
   ;; if set as a keyword, apply to all transits according to the keyword
   fitdilute = lonarr(ntran>1) + keyword_set(fitdilute)
endif else if n_elements(fitdilute) eq ntran then begin
   ;; do nothing, user has already fully specified FITDILUTE for each light curve
endif else begin
   ;; this is an error
   printandlog, 'FITDILUTE (' + strtrim(n_elements(fitdilute),2) + ') must be an NTRAN (' + strtrim(ntran,2) + ') array', logname
   return, -1
endelse

if ntran eq 0 then begin
   exptime = [1]
   ninterp = [1]
endif else begin
   if n_elements(ninterp) eq 1 and n_elements(exptime) eq 1 then begin
      exptime = dblarr(ntran) + exptime
      ninterp = dblarr(ntran) + ninterp
   endif else if n_elements(ninterp) eq ntran and n_elements(exptime) eq ntran then begin
      ;; do nothing (use exptime and ninterp as is)
   endif else if n_elements(ninterp) eq 0 and n_elements(exptime) eq 0 then begin
      if n_elements(longcadence) eq 0 then begin
         exptime = dblarr(ntran) + 1
         ninterp = dblarr(ntran) + 1
      endif else if n_elements(longcadence) eq 1 then begin
         if longcadence[0] eq 1 then begin
            exptime = dblarr(ntran) + 29.425d0
            ninterp = dblarr(ntran) + 10
         endif else begin
            exptime = dblarr(ntran) + 1
            ninterp = dblarr(ntran) + 1
         endelse
      endif else if n_elements(longcadence) ne ntran then begin
         printandlog, 'LONGCADENCE must be byte or an NTRANSITS (' + strtrim(ntran,2) + ') byte array', logname
         return, -1
      endif else begin
         exptime = dblarr(ntran) + 29.425d0
         ninterp = dblarr(ntran) + 1
         match = where(longcadence)
         if match[0] ne -1 then ninterp[match] = 10
      endelse
   endif else begin
      printandlog, 'NINTERP and EXPTIME must be unspecified or an NTRANSITS (' + strtrim(ntran,2) + ') array', logname
      return, -1
   endelse
endelse

bad = where(ninterp le 0 or exptime le 0, nbad)
if nbad ne 0 then begin
   printandlog, 'NINTERP and EXPTIME must be greater than zero',logname
   return, -1
endif

if nplanets ge 1 and ntran ge 1 then begin
   ;; some error checking on TTVs
   if n_elements(ttvs) eq 0 then ttvs = bytarr(ntran,nplanets) $
   else if n_elements(ttvs) eq 1 then ttvs = bytarr(ntran,nplanets)+ttvs[0] $
   else if n_elements(ttvs) ne ntran*nplanets then begin
      printandlog, 'TTVs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      return, -1
   endif
   if nplanets gt 1 then begin
      if max(total(ttvs,2)) gt 1 then begin
         printandlog, 'TTVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
   
   ;; some error checking on TIVs
   if n_elements(tivs) eq 0 then tivs = bytarr(ntran,nplanets) $
   else if n_elements(tivs) eq 1 then tivs = bytarr(ntran,nplanets)+tivs[0] $
   else if n_elements(tivs) ne ntran*nplanets then begin
      printandlog, 'TIVs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      return, -1
   endif
   if nplanets gt 1 then begin
      if max(total(tivs,2)) gt 1 then begin
         printandlog, 'TIVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
   
   ;; some error checking on TDELTAVs
   if n_elements(tdeltavs) eq 0 then tdeltavs = bytarr(ntran,nplanets) $
   else if n_elements(tdeltavs) eq 1 then tdeltavs = bytarr(ntran,nplanets)+tdeltavs[0] $
   else if n_elements(tdeltavs) ne ntran*nplanets then begin
      printandlog, 'Tdeltavs must be an NTRANSITSxNPLANETS (' + string(ntran,nplanets,format='(i,"x",i)') + ') array', logname
      return, -1
   endif
   if nplanets gt 1 then begin
      if max(total(tdeltavs,2)) gt 1 then begin
         printandlog, 'TDELTAVs applied to more than one planet. This is unlikely to be what you want.', logname
      endif
   endif
endif else begin
   ttvs = 0B
   tivs = 0B
   tdeltavs = 0B
endelse

if n_elements(rejectflatmodel) ne ntran and n_elements(rejectflatmodel) ne 0 and ntran gt 0 then begin
   printandlog, 'REJECTFLATMODEL must be an NTRANSITS element array', logname
   return, -1
end
if n_elements(rejectflatmodel) eq 0 then begin
   if ntran gt 0 then rejectflatmodel = bytarr(ntran) $
   else rejectflatmodel = [0B]
endif

if n_elements(fitspline) ne ntran and n_elements(fitspline) ne 0 and ntran gt 0 then begin
   printandlog, 'FITSPLINE has ' + strtrim(n_elements(fitspline),2) + ' elements; must be an NTRANSITS (' + strtrim(ntran,2) + ') element array', logname
   return, -1
end
if n_elements(fitspline) eq 0 then begin
   if ntran gt 0 then fitspline = bytarr(ntran) $
   else fitspline = [0B]
endif

if n_elements(splinespace) ne ntran and n_elements(splinespace) ne 0 and ntran gt 0 then begin
   printandlog, 'SPLINESPACE must be an NTRANSITS element array', logname
   return, -1
end
if n_elements(splinespace) eq 0 then begin
   if ntran gt 0 then splinespace = dblarr(ntran) + 0.75d0 $
   else splinespace = [0.75d0]
endif

if n_elements(fitwavelet) ne ntran and n_elements(fitwavelet) ne 0 and ntran gt 0 then begin
   printandlog, 'FITWAVELET must be an NTRANSITS element array', logname
   return, -1
end
if n_elements(fitwavelet) eq 0 then begin
   if ntran gt 0 then fitwavelet = bytarr(ntran) $
   else fitwavelet = [0B]
endif

;; if RVPATH empty, don't use any telescope
if rvpath eq '' then begin
   ntel = 0
   rvpath = ''
   if keyword_set(fitslope) or keyword_set(fitquad) then begin
      printandlog, "WARNING: RV path not specified. Ignoring FITSLOPE and FITQUAD keywords." , logname
   endif
endif else begin
   rvfiles = file_search(rvpath,count=ntel)
   ;; was it specifed and not found?
   if ntel eq 0 then begin
      printandlog, "RV path (" + rvpath + ") not found! Make sure the file exists or remove the argument to proceed without it." , logname
      return, -1
   endif
endelse

;; same for DT path
if dtpath ne '' then begin
   dtfiles = file_search(dtpath,count=ndt)
   if ndt eq 0 then message, 'DTPATH specified (' + dtpath + ') but no files found'
endif else begin
   dtfiles = ['']
   ndt = 0
endelse

if n_elements(circular) ne nplanets or nplanets eq 0 then circular = bytarr(nplanets>1)
if n_elements(rossiter) ne nplanets or nplanets eq 0 then rossiter = bytarr(nplanets>1)
if n_elements(novcve) ne nplanets or nplanets eq 0 then novcve = bytarr(nplanets>1)
if n_elements(nochord) ne nplanets or nplanets eq 0 then nochord = bytarr(nplanets>1)
if n_elements(fitsign) ne nplanets or nplanets eq 0 then fitsign = bytarr(nplanets>1)+1B
if n_elements(fitbeam) ne nplanets or nplanets eq 0 then fitbeam = bytarr(nplanets>1)
if n_elements(derivebeam) ne nplanets or nplanets eq 0 then derivebeam = bytarr(nplanets>1)

if n_elements(nvalues) ne 0 then value = dblarr(nvalues) $
else value = 0d0
nsteps = n_elements(value)

if n_elements(fittran) eq 0 then begin
   if ntran eq 0 then fittran = bytarr(nplanets>1) $
   else fittran = bytarr(nplanets>1)+1B
endif

if n_elements(fitrv) eq 0 then begin
   if ntel eq 0 then fitrv = bytarr(nplanets>1) $
   else fitrv = bytarr(nplanets>1)+1B
endif

if n_elements(fitdt) eq 0 then begin
   if ndt eq 0 then fitdt = bytarr(nplanets>1) $
   else fitdt = bytarr(nplanets>1)+1B
endif

if n_elements(fittt) eq 0 then begin
   fittt = bytarr(nplanets>1)
endif

if n_elements(fitlogmp) eq 0 then fitlogmp = bytarr(nplanets>1)

if n_elements(chen) ne nplanets or nplanets eq 0 then chen = fittran xor fitrv
if n_elements(i180) ne nplanets or nplanets eq 0 then i180 = bytarr(nplanets>1)

if n_elements(fittran) ne nplanets and nplanets ne 0 then begin
   printandlog, "FITTRAN must have NPLANETS elements", logname
   return, -1
endif

if n_elements(fitrv) ne nplanets and nplanets ne 0 then begin
   printandlog, "FITRV must have NPLANETS elements", logname
   return, -1
endif

if nplanets ne 0 then begin
   if (where((~fitrv) and (~fittran)))[0] ne -1 and astrompath eq '' then begin
      printandlog, 'Either a transit or RV must be fit for each planet', logname
      return, -1
   end
endif

;; was rvpath specified but no planets are fit? (ignore it)
junk = where(fitrv,nrvfit)
if rvpath ne '' and nrvfit eq 0 then begin
   if ~keyword_set(silent) then printandlog, 'WARNING: an RVPATH was specified, but no RV planets are fit. Ignoring the supplied RVs'
   if ~keyword_set(silent) then printandlog, 'Remove RVPATH or set at least one planet in FITRV to 1 to remove this message.'
   ntel = 0
endif

if rvpath eq '' and nrvfit ne 0 then begin
   printandlog, 'ERROR: a RVPATH was not specified, but RVs were requested to be fit (FITRV != 0).'
   printandlog, 'Remove FITRV or specify a RVPATH'
   return, -1
endif

;; was tranpath specified but no planets are fit? (ignore it)
junk = where(fittran,ntranfit)
if tranpath ne '' and ntranfit eq 0 then begin
   if ~keyword_set(silent) then printandlog, 'WARNING: a TRANPATH was specified, but no transits are fit. Ignoring the supplied transits'
   if ~keyword_set(silent) then printandlog, 'Remove TRANPATH or set at least one planet in FITTRAN to 1 to remove this message.'
   ntran = 0
endif

if tranpath eq '' and ntranfit ne 0 then begin
   printandlog, 'ERROR: a TRANPATH was not specified, but transits were requested to be fit (FITTRAN != 0).'
   printandlog, 'Remove FITTRAN or specify a TRANPATH'
   return, -1
endif

;; each parameter is a structure, as defined here
parameter = create_struct('value',value,$     ;; its numerical value
                          'prior',0d0,$       ;; its prior value
                          'priorwidth',!values.d_infinity,$ ;; its prior width (infinity => no constraint)
                          'lowerbound',-!values.d_infinity,$ ;; values lower than this have zero likelihood
                          'upperbound',!values.d_infinity,$ ;; values higher than this have zero likelihood
                          'label','',$        ;; what do I call it?
                          'cgs',1d0,$         ;; multiply value by this to convert to cgs units
                          'scale',0d0,$       ;; scale for amoeba
                          'prior_new',ptr_new(/allocate_heap),$
                          'latex','',$        ;; latex label for the latex table
                          'userchanged',0B,$  ;; the user supplied a prior that impacts this parameter
                          'description','',$  ;; a verbose description
                          'unit','',$         ;; units ("Radians" triggers special angular behavior)
                          'medvalue','',$     ;; median value
                          'upper','',$        ;; upper error bar (68% confidence)
                          'lower','',$        ;; lower error bar (68% confidence)
                          'scinote','',$      ;; text for scientific notation if applicable
                          'best',!values.d_nan,$ ;; best-fit parameter
                          'fit', 0B,$         ;; If true, step in this parameter during the MCMC
                          'derive',1B)        ;; If true, quote this parameter in the final table

null_parameter = create_struct('value',0L,$                  ;; its numerical value
                               'prior',0d0,$                 ;; its prior value
                               'priorwidth',!values.d_infinity,$ ;; its prior width (infinity => no constraint)
                               'lowerbound',-!values.d_infinity,$ ;; values lower than this have zero likelihood
                               'upperbound',!values.d_infinity,$  ;; values higher than this have zero likelihood
                               'label','',$                       ;; what do I call it?
                               'cgs',1d0,$                        ;; multiply value by this to convert to cgs units
                               'scale',0d0,$                      ;; scale for amoeba
                               'prior_new',ptr_new(/allocate_heap),$
                               'latex','',$        ;; latex label for the latex table
                               'userchanged',0B,$  ;; the user supplied a prior that impacts this parameter
                               'description','',$  ;; a verbose description
                               'unit','',$         ;; units ("Radians" triggers special angular behavior)
                               'medvalue','',$     ;; median value
                               'upper','',$        ;; upper error bar (68% confidence)
                               'lower','',$        ;; lower error bar (68% confidence)
                               'best',!values.d_nan,$ ;; best-fit parameter
                               'fit', 0B,$            ;; If true, step in this parameter during the MCMC
                               'derive',1B)           ;; If true, quote this parameter in the final table

                          
;; define each parameter (derived, fitted, and ignored)
logmstar = parameter
logmstar.value = 0d0
logmstar.unit = ''
logmstar.description = 'log Mass'
logmstar.latex = '\log{M_*/M_{\sun}}'
logmstar.label = 'logmstar'
logmstar.cgs = !values.d_nan
logmstar.fit = 1
logmstar.derive=0
logmstar.scale = 1

absks = parameter
absks.value = 10d0
absks.unit = 'mag'
absks.description = 'Absolute Ks-band mag'
absks.latex = 'K_S'
absks.label = 'absks'
absks.cgs = !values.d_nan
absks.fit = 0
absks.derive=0
absks.scale = 0.1

appks = parameter
appks.value = 10d0
appks.unit = 'mag'
appks.description = 'Apparent Ks-band mag'
appks.latex = 'k_S'
appks.label = 'appks'
appks.cgs = !values.d_nan
appks.fit = 0
appks.derive=0
appks.scale = 0.1

mstar = parameter
mstar.value = 10^logmstar.value
mstar.unit = '\msun'
mstar.description = 'Mass'
mstar.latex = 'M_*'
mstar.label = 'mstar'
mstar.cgs = 1.9891d33

logg = parameter
logg.value = 4.44d0
logg.unit = 'cgs'
logg.description = 'Surface gravity'
logg.latex = '\log{g}'
logg.label = 'logg'
logg.scale = 0.3d0

teff = parameter
teff.value = 5778d0
teff.unit = 'K'
teff.description = 'Effective temperature'
teff.latex = 'T_{\rm eff}'
teff.label = 'teff'
teff.fit=1
teff.scale = 500d0

teffsed = parameter
teffsed.value = 5778d0
teffsed.unit = 'K'
teffsed.description = 'Effective temperature'
teffsed.latex = 'T_{\rm eff,SED}'
teffsed.label = 'teffsed'
teffsed.fit=0
teffsed.derive=0
teffsed.scale = 500d0

feh = parameter
feh.value = 0d0
feh.description = 'Metallicity'
feh.latex = '[{\rm Fe/H}]'
feh.label = 'feh'
feh.unit = 'dex'
feh.fit=1
feh.scale = 0.5d0

fehsed = parameter
fehsed.value = 0d0
fehsed.description = 'Metallicity'
fehsed.latex = '[{\rm Fe/H}]_{SED}'
fehsed.label = 'fehsed'
fehsed.unit = 'dex'
fehsed.fit=0
fehsed.derive=0
fehsed.scale = 0.5d0

initfeh = parameter
initfeh.value = 0d0
initfeh.description = 'Initial Metallicity'
initfeh.latex = '[{\rm Fe/H}]_{0}'
initfeh.label = 'initfeh'
initfeh.scale = 0.5d0

Av = parameter
Av.description = 'V-band extinction'
Av.latex = 'A_V'
Av.label = 'Av'
Av.unit = 'mag'
Av.fit = 0
Av.derive = 0
Av.scale = 0.3d0

alpha = parameter
alpha.description = 'Alpha abundance'
alpha.latex = '\alpha'
alpha.label = 'alpha'
alpha.fit = 0
alpha.derive = 0
alpha.scale = 0.3d0

;Ma = parameter
;Ma.unit = ''
;Ma.description = 'Apparent V-band Magnitude'
;Ma.latex = 'M_a'
;Ma.label = 'Ma'
;Ma.fit=1

;Mv = parameter
;Mv.unit = ''
;Mv.description = 'Absolute V-band Magnitude'
;Mv.latex = 'M_v'
;Mv.label = 'Mv'

distance = parameter
distance.unit = 'pc'
distance.description = 'Distance'
distance.latex = 'd'
distance.label = 'distance'
distance.cgs = 3.08567758d18 ;; cm/pc
distance.scale = 100
distance.value = 10
distance.fit = 0B
distance.derive = 0B
if nastrom gt 0 then begin
   distance.fit = 1d0
   distance.derive = 1B
endif

fbol = parameter
fbol.unit = 'cgs'
fbol.description = 'Bolometric Flux'
fbol.latex = 'F_{Bol}'
fbol.label = 'fbol'
fbol.cgs = 1d0
fbol.derive = 0d0

ra = parameter
ra.value = 0d0
ra.scale = 0.1d0
ra.unit = 'ICRS deg'
ra.description = 'Right Ascension'
ra.latex = '\alpha_{ICRS}'
ra.label = 'ra'
ra.cgs = 3600d0*180d0/!dpi ;; rad/deg
ra.derive = 0B

dec = parameter
dec.value = 0d0
dec.scale = 0.1d0
dec.unit = 'ICRS deg'
dec.description = 'Declination'
dec.latex = '\delta_{ICRS}'
dec.label = 'dec'
dec.cgs = 3600d0*180d0/!dpi ;; rad/deg
dec.derive = 0B

pmra = parameter
pmra.value = 0d0
pmra.scale = 1000d0
pmra.unit = 'mas/yr'
pmra.description = 'RA Proper Motion'
pmra.latex = '\mu_{\alpha*}'
pmra.label = 'pmra'
pmra.cgs = 3600d3*180d0/!dpi*86400d0*365.25d0 ;; rad/deg
pmra.derive = 0B

pmdec = parameter
pmdec.value = 0d0
pmdec.scale = 1000d0
pmdec.unit = 'mas/yr'
pmdec.description = 'Dec Proper Motion'
pmdec.latex = '\mu_{\delta}'
pmdec.label = 'pmdec'
pmdec.cgs = 3600d3*180d0/!dpi*86400d0*365.25d0 ;; rad/deg
pmdec.derive = 0B

rvabs = parameter
rvabs.value = 0d0
rvabs.scale = 10000d0
rvabs.unit = 'm/s'
rvabs.description = 'Absolute RV'
rvabs.latex = '\gamma_{\rm abs}'
rvabs.label = 'rvabs'
rvabs.cgs = 100d0 ;; m/s
rvabs.derive = 0B

astromscale = parameter
astromscale.value = 1d0
astromscale.scale = 10d0
astromscale.unit = ''
astromscale.description = 'Astrometric Error Scaling'
astromscale.latex = '\sigma_{Astrom}'
astromscale.label = 'astromscale'
astromscale.cgs = 0
if nastrom gt 0 then begin
   astromscale.fit = 1B
   astromscale.derive=1B
endif else astromscale.derive = 0B

raoffset = parameter
raoffset.value = 0d0
raoffset.scale = 1d-2
raoffset.unit = 'degrees'
raoffset.description = 'RA offset'
raoffset.latex = '\Delta\alpha'
raoffset.label = 'raoffset'
raoffset.cgs = 0 ;;
raoffset.derive = 0B

decoffset = parameter
decoffset.value = 0d0
decoffset.scale = 1d-2
decoffset.unit = 'degrees'
decoffset.description = 'Dec offset'
decoffset.latex = '\Delta\delta'
decoffset.label = 'decoffset'
decoffset.cgs = 1d0 ;;
decoffset.derive = 0B

phottobary = parameter
phottobary.value = 1d0
phottobary.scale = 10d0
phottobary.unit = ''
phottobary.description = 'Photocenter/Barycenter'
phottobary.latex = 'p/\rho'
phottobary.label = 'phottobary'
phottobary.cgs = 1d0
phottobary.derive = 0B

parallax = parameter
parallax.value = 0d0
parallax.scale = 100d0
parallax.unit = 'mas'
parallax.description = 'Parallax'
parallax.latex = '\varpi'
parallax.label = 'parallax'
parallax.cgs = 3600d3*180d0/!dpi ;; rad/mas
parallax.derive=0B
parallax.fit = 0B
if nastrom gt 0 then begin
   parallax.derive=1B
   fbol.derive = 1B
endif

gamma = parameter
gamma.unit = 'm/s'
gamma.description = 'Relative RV Offset'
gamma.latex = '\gamma_{\rm rel}'
gamma.label = 'gamma'
gamma.cgs = 100d0
if ntel gt 0 and nplanets gt 0 then gamma.fit=1 $
else gamma.derive=0
gamma.scale=5000d0

slope = parameter
slope.unit = 'm/s/day'
slope.description = 'RV slope'
slope.latex = '\dot{\gamma}'
slope.label = 'slope'
slope.cgs = 100d0/86400d0
if (keyword_set(fitslope) or keyword_set(fitquad)) and (rvpath ne '') then slope.fit = 1 $
else slope.derive=0
slope.scale=1d0

quad = parameter
quad.unit = 'm/s/day^2'
quad.description = 'RV quadratic term'
quad.latex = '\ddot{\gamma}'
quad.label = 'quad'
quad.cgs = 100d0/86400d0^2
if keyword_set(fitquad) and rvpath ne '' then quad.fit=1 $
else quad.derive=0
quad.scale = 1d0

rstar = parameter
rstar.value = 1d0
rstar.unit = '\rsun'
rstar.description = 'Radius'
rstar.latex = 'R_*'
rstar.label = 'rstar'
rstar.cgs = 6.955d10
rstar.fit = 1
rstar.scale = 0.5d0

rstarsed = parameter
rstarsed.value = 1d0
rstarsed.unit = '\rsun'
rstarsed.description = 'Radius'
rstarsed.latex = 'R_{*,SED}'
rstarsed.label = 'rstarsed'
rstarsed.cgs = 6.955d10
rstarsed.fit = 0
rstarsed.derive = 0
rstarsed.scale = 0.5d0

age = parameter
age.value = 4.603d0
age.unit = 'Gyr'
age.description = 'Age'
age.latex = 'Age'
age.label = 'age'
age.cgs = 3600d0*24d0*365.242d0*1d9
age.derive = 0
age.scale = 3d0

eep = parameter
eep.value = 398.668d0 ;; solar EEP (PARSEC)
eep.value = 354.1661d0 ;; solar EEP (MIST)
eep.unit = ''
eep.description = 'Equal Evolutionary Phase'
eep.latex = 'EEP'
eep.label = 'eep'
eep.scale = 50d0

lstar = parameter
lstar.unit = '\lsun'
lstar.description = 'Luminosity'
lstar.latex = 'L_*'
lstar.label = 'lstar'
lstar.cgs = 3600d0*24d0*365.242d0*1d9

rhostar = parameter
rhostar.unit = 'cgs'
rhostar.description = 'Density'
rhostar.latex = '\rho_*'
rhostar.label = 'rhostar'

tc = parameter
tc.unit = '\bjdtdb'
tc.description = 'Time of conjunction'
tc.latex = 'T_C'
tc.label = 'tc'
tc.cgs = 86400d0
if nplanets eq 0 then tc.derive=0 $
else tc.fit = 1
tc.scale = 0.1

tt = parameter
tt.unit = '\bjdtdb'
tt.description = 'Time of min proj sep'
tt.latex = 'T_T'
tt.label = 'tt'
tt.cgs = 86400d0
if ntran eq 0 then tt.derive=0
tt.scale = 0.1
tt.derive = 0B

t0 = parameter
t0.unit = '\bjdtdb'
t0.description = 'Optimal conj time'
t0.latex = 'T_0'
t0.label = 't0'
t0.cgs = 86400d0
if nplanets eq 0 then t0.derive=0

logp = parameter
logp.unit = ''
logp.description = 'Log of Period'
logp.latex = '\log{P}'
logp.label = 'logp'
logp.cgs = !values.d_nan
if nplanets ne 0 then logp.fit = 1
logp.derive=0
logp.scale = 0.01d0

qesinw = parameter
qesinw.unit = ''
qesinw.description = ''
qesinw.latex = 'e^{1/4}\sin{\omega_*}'
qesinw.label = 'qesinw'
qesinw.cgs = !values.d_nan
qesinw.scale = 0.1d0
qesinw.derive = 0

qecosw = parameter
qecosw.unit = ''
qecosw.description = ''
qecosw.latex = 'e^{1/4}\cos{\omega_*}'
qecosw.label = 'qecosw'
qecosw.cgs = !values.d_nan
qecosw.scale = 0.1d0
qecosw.derive = 0

sesinw = parameter
sesinw.unit = ''
sesinw.description = ''
sesinw.latex = '\sqrt{e}\sin{\omega_*}'
sesinw.label = 'sesinw'
sesinw.cgs = !values.d_nan
sesinw.scale = 0.1d0
sesinw.derive = 0

secosw = parameter
secosw.unit = ''
secosw.description = ''
secosw.latex = '\sqrt{e}\cos{\omega_*}'
secosw.label = 'secosw'
secosw.cgs = !values.d_nan
secosw.scale = 0.1d0
secosw.derive = 0

vcve = parameter
vcve.unit = ''
vcve.description = 'Scaled velocity'
vcve.latex = 'V_c/V_e'
vcve.label = 'vcve'
vcve.cgs = !values.d_nan
vcve.scale = 1d0
vcve.value = 0.99d0
if nplanets eq 0 then vcve.derive=0B

chord = parameter
chord.unit = ''
chord.description = 'Transit chord'
chord.latex = '((1-R_P/R_*)^2-b^2)^{1/2}'
chord.label = 'chord'
chord.cgs = 1d0
chord.scale = 1d0
chord.value = 1d0
chord.derive = 0

sign = parameter
sign.unit = ''
sign.description = ''
sign.latex = 'sign'
sign.label = 'sign'
sign.cgs = !values.d_nan
sign.scale = 1
sign.value = 1.5
sign.derive = 0

sign2 = parameter
sign2.unit = ''
sign2.description = ''
sign2.latex = 'sign2'
sign2.label = 'sign2'
sign2.cgs = !values.d_nan
sign2.scale = 1
sign2.value = 1.5
sign2.derive = 0

if keyword_set(eprior4) then begin
   ;; step in qesinw (prior favoring lower e)
   qesinw.fit = 1
   qecosw.fit = 1
endif else begin
   ;; step in sqrt(e)sinw (uniform prior)
   if nplanets ne 0 then begin
      sesinw.fit = 1
      secosw.fit = 1
   endif
endelse

esinw = parameter
esinw.unit = ''
esinw.description = ''
esinw.latex = 'e\sin{\omega_*}'
esinw.label = 'esinw'
esinw.cgs = !values.d_nan
if nplanets eq 0 then esinw.derive=0

ecosw = parameter
ecosw.unit = ''
ecosw.description = ''
ecosw.latex = 'e\cos{\omega_*}'
ecosw.label = 'ecosw'
ecosw.cgs = !values.d_nan
if nplanets eq 0 then ecosw.derive=0

e = parameter
e.unit = ''
e.description = 'Eccentricity'
e.latex = 'e'
e.label = 'e'
e.cgs = !values.d_nan
if nplanets eq 0 then e.derive=0

omega = parameter
omega.unit = 'Radians'
omega.description = 'Arg of periastron'
omega.latex = '\omega_*'
omega.label = 'omega'
omega.cgs = 1d0
omega.derive = 0

lsinw = parameter
lsinw.unit = ''
lsinw.description = 'L*Sine of Argument of Periastron'
lsinw.latex = 'L\sin\{\omega_*}'
lsinw.label = 'lsinw'
lsinw.cgs = 1d0
lsinw.value = 0d0
lsinw.scale = 1d0
lsinw.derive = 0

lsinw2 = parameter
lsinw2.unit = ''
lsinw2.description = 'L*Sine of Argument of Periastron'
lsinw2.latex = 'sign(vcve-1)*L\sin\{\omega_*}'
lsinw2.label = 'lsinw2'
lsinw2.cgs = 1d0
lsinw2.value = 0d0
lsinw2.scale = 1d0
lsinw2.derive = 0

lcosw = parameter
lcosw.unit = ''
lcosw.description = 'L*Cos of arg of periastron'
lcosw.latex = 'L\cos{\omega_*}'
lcosw.label = 'lcosw'
lcosw.cgs = 1d0
lcosw.value = 0d0
lcosw.scale = 1d0
lcosw.derive = 0

omegadeg = parameter
omegadeg.unit = 'Degrees'
omegadeg.description = 'Arg of periastron'
omegadeg.latex = '\omega_*'
omegadeg.label = 'omegadeg'
omegadeg.cgs = !dpi/180d0
if nplanets eq 0 then omegadeg.derive=0

bigomega = parameter
bigomega.value = 0d0
bigomega.scale = !dpi
bigomega.unit = 'Radians'
bigomega.description = 'Long of asc node'
bigomega.latex = '\Omega'
bigomega.label = 'bigomega'
bigomega.cgs = 1d0
bigomega.derive = 0B

lsinbigomega = parameter
lsinbigomega.value = 0d0
lsinbigomega.scale = 0.5
lsinbigomega.unit = ''
lsinbigomega.description = 'l*Sine of Longitude of ascending node'
lsinbigomega.latex = 'l\sin{\Omega}'
lsinbigomega.label = 'lsinbigomega'
lsinbigomega.cgs = 1d0
lsinbigomega.derive = 0B

lcosbigomega = parameter
lcosbigomega.value = 0d0
lcosbigomega.scale = 0.5
lcosbigomega.unit = ''
lcosbigomega.description = 'l*Cosine of Longitude of ascending node'
lcosbigomega.latex = 'l\cos{\Omega}'
lcosbigomega.label = 'lcosbigomega'
lcosbigomega.cgs = 1d0
lcosbigomega.derive = 0B

bigomegadeg = parameter
bigomegadeg.value = 0d0
bigomegadeg.scale = 180d0
bigomegadeg.unit = 'Degrees'
bigomegadeg.description = 'Longitude of ascending node'
bigomegadeg.latex = '\Omega'
bigomegadeg.label = 'bigomegadeg'
bigomegadeg.cgs = !dpi/180d0
bigomegadeg.derive=0B
bigomegadeg.fit=0B
if nastrom gt 0 then begin
   lsinbigomega.fit=1B
   lcosbigomega.fit=1B
   bigomegadeg.derive=1B
endif else begin
   lsinbigomega.fit=0B
   lcosbigomega.fit=0B
endelse
   
p = parameter
p.value = 0.1d0
p.unit = ''
p.description = 'Radius of planet in stellar radii'
p.latex = 'R_P/R_*'
p.label = 'p'
p.scale = 1d-1
if nplanets eq 0 then p.derive = 0 $
else p.fit = 1

ar = parameter
ar.unit = ''
ar.description = 'Semi-major axis in stellar radii'
ar.latex = 'a/R_*'
ar.label = 'ar'
if nplanets eq 0 then ar.derive=0

cosi = parameter
cosi.value = 0d0
cosi.unit = ''
cosi.description = 'Cos of inclination'
cosi.latex = '\cos{i}'
cosi.label = 'cosi'
if nplanets ne 0 then cosi.fit = 1
cosi.scale = 0.1
cosi.derive = 0

inc = parameter
inc.unit = 'Radians'
inc.description = 'Inclination'
inc.latex = 'i'
inc.label = 'i'
inc.derive = 0

ideg = parameter
ideg.unit = 'Degrees'
ideg.description = 'Inclination'
ideg.latex = 'i'
ideg.label = 'ideg'
ideg.cgs = !dpi/180d0
if nplanets eq 0 then ideg.derive = 0

tp = parameter
tp.unit = '\bjdtdb'
tp.description = 'Time of Periastron'
tp.latex = 'T_P'
tp.label = 'tp'
tp.cgs = 86400d0
if nplanets eq 0 then tp.derive = 0

td = parameter
td.unit = '\bjdtdb'
td.description = 'Time of desc node'
td.latex = 'T_D'
td.label = 'td'
td.cgs = 86400d0
if nplanets eq 0 then td.derive = 0

ta = parameter
ta.unit = '\bjdtdb'
ta.description = 'Time of asc node'
ta.latex = 'T_A'
ta.label = 'ta'
ta.cgs = 86400d0
if nplanets eq 0 then ta.derive = 0

ts = parameter
ts.unit = '\bjdtdb'
ts.description = 'Time of eclipse'
ts.latex = 'T_S'
ts.label = 'ts'
ts.cgs = 86400d0
if nplanets eq 0 then ts.derive = 0

phase = parameter
phase.unit = ''
phase.description = 'Phase of inferior conjunction (primary transit)'
phase.latex = '\phi'
phase.label = 'phase'
phase.derive = 0

tfwhm = parameter
tfwhm.unit = 'days'
tfwhm.description = 'FWHM transit duration'
tfwhm.latex = 'T_{FWHM}'
tfwhm.label = 'tfwhm'
tfwhm.cgs = 86400d0
if nplanets eq 0 then tfwhm.derive = 0

t14 = parameter
t14.unit = 'days'
t14.description = 'Total transit duration'
t14.latex = 'T_{14}'
t14.label = 't14'
t14.cgs = 86400d0
if nplanets eq 0 then t14.derive = 0

tau = parameter
tau.unit = 'days'
tau.description = 'In/egress transit duration'
tau.latex = '\tau'
tau.label = 'tau'
tau.cgs = 86400d0
if nplanets eq 0 then tau.derive = 0

tfwhms = parameter
tfwhms.unit = 'days'
tfwhms.description = 'FWHM eclipse duration'
tfwhms.latex = 'T_{S,FWHM}'
tfwhms.label = 'tfwhms'
tfwhms.cgs = 86400d0
if nplanets eq 0 then tfwhms.derive = 0

t14s = parameter
t14s.unit = 'days'
t14s.description = 'Total eclipse duration'
t14s.latex = 'T_{S,14}'
t14s.label = 't14s'
t14s.cgs = 86400d0
if nplanets eq 0 then t14s.derive = 0

taus = parameter
taus.unit = 'days'
taus.description = 'In/egress eclipse duration'
taus.latex = '\tau_S'
taus.label = 'taus'
taus.cgs = 86400d0
if nplanets eq 0 then taus.derive = 0

eclipsedepth = parameter
eclipsedepth.unit = 'ppm'
eclipsedepth.description = 'Measured eclipse depth'
eclipsedepth.latex = '\delta_{S}'
eclipsedepth.label = 'eclipsedepth'
eclipsedepth.cgs = 1d6
eclipsedepth.derive = 0

;eclipsedepth36 = parameter
;eclipsedepth36.unit = 'ppm'
;eclipsedepth36.description = 'Blackbody eclipse depth at 3.6$\mu$m'
;eclipsedepth36.latex = '\delta_{S,3.6\mu m}'
;eclipsedepth36.label = 'eclipsedepth36'
;eclipsedepth36.cgs = 1d6
;if nplanets eq 0 then eclipsedepth36.derive = 0

;eclipsedepth45 = parameter
;eclipsedepth45.unit = 'ppm'
;eclipsedepth45.description = 'Blackbody eclipse depth at 4.5$\mu$m'
;eclipsedepth45.latex = '\delta_{S,4.5\mu m}'
;eclipsedepth45.label = 'eclipsedepth45'
;eclipsedepth45.cgs = 1d6
;if nplanets eq 0 then eclipsedepth45.derive = 0

eclipsedepth25 = parameter
eclipsedepth25.unit = 'ppm'
eclipsedepth25.description = 'BB eclipse depth at 2.5$\mu$m'
eclipsedepth25.latex = '\delta_{S,2.5\mu m}'
eclipsedepth25.label = 'eclipsedepth25'
eclipsedepth25.cgs = 1d6
if nplanets eq 0 then eclipsedepth25.derive = 0

eclipsedepth50 = parameter
eclipsedepth50.unit = 'ppm'
eclipsedepth50.description = 'BB eclipse depth at 5.0$\mu$m'
eclipsedepth50.latex = '\delta_{S,5.0\mu m}'
eclipsedepth50.label = 'eclipsedepth50'
eclipsedepth50.cgs = 1d6
if nplanets eq 0 then eclipsedepth50.derive = 0

eclipsedepth75 = parameter
eclipsedepth75.unit = 'ppm'
eclipsedepth75.description = 'BB eclipse depth at 7.5$\mu$m'
eclipsedepth75.latex = '\delta_{S,7.5\mu m}'
eclipsedepth75.label = 'eclipsedepth75'
eclipsedepth75.cgs = 1d6
if nplanets eq 0 then eclipsedepth75.derive = 0

delta = parameter
;delta.unit = 'frac'
delta.description = '$\left(R_P/R_*\right)^2$'
delta.latex = '\delta'
delta.label = 'delta'
if nplanets eq 0 then delta.derive = 0

depth = parameter
depth.description = 'Flux decrement at mid transit'
depth.latex = 'Depth'
depth.label = 'depth'
if nplanets eq 0 then depth.derive = 0

dr = parameter
dr.description = 'Separation at mid transit'
dr.latex = 'd/R_*'
dr.label = 'dr'
if nplanets eq 0 then dr.derive = 0

lambda = parameter
lambda.unit = 'Radians'
lambda.description = 'Projected Spin-orbit alignment'
lambda.latex = '\lambda'
lambda.label = 'lambda'
lambda.scale = !dpi
lambda.derive = 0

lsinlambda = parameter
lsinlambda.value = 0d0
lsinlambda.unit = ''
lsinlambda.description = 'lSine of Projected Spin-orbit alignment'
lsinlambda.latex = 'l\sin{\lambda}'
lsinlambda.label = 'lsinlambda'
lsinlambda.scale = 0.5
lsinlambda.derive = 0

lcoslambda = parameter
lcoslambda.value = 0d0
lcoslambda.unit = ''
lcoslambda.description = 'l*Cosine of Projected Spin-orbit alignment'
lcoslambda.latex = 'l\cos{\lambda}'
lcoslambda.label = 'lcoslambda'
lcoslambda.scale = 0.5
lcoslambda.derive = 0

lambdadeg = parameter
lambdadeg.unit = 'Degrees'
lambdadeg.description = 'Projected Spin-orbit alignment'
lambdadeg.latex = '\lambda'
lambdadeg.label = 'lambdadeg'
lambdadeg.cgs = !dpi/180d0
lambdadeg.derive = 0

vsini = parameter
vsini.unit = 'm/s'
vsini.description = 'Projected rotational velocity'
vsini.latex = 'vsinI_*'
vsini.label = 'vsini'
vsini.cgs = 100d0
vsini.derive = 0
vsini.scale = 5d3

vline = parameter
vline.unit = 'm/s'
vline.description = 'Unbroadened line width'
vline.latex = 'V_{\rm line}'
vline.label = 'vline'
vline.cgs = 1000d0
vline.derive = 0
vline.value = 4d3
vline.scale = 2d3

dtscale = parameter
dtscale.description = 'Doppler Tomography Error scaling'
dtscale.latex = '\sigma_{DT}'
dtscale.label = 'dtscale'
dtscale.derive = 0
dtscale.scale = 1d2
dtscale.value = 1d0

logmp = parameter
logmp.description = 'Log of mp in solar units'
logmp.latex = '\log{M_P}'
logmp.label = 'logmp'
logmp.derive=0

k = parameter
k.value = 10d0
k.unit = 'm/s'
k.description = 'RV semi-amplitude'
k.latex = 'K'
k.label = 'k'
k.cgs = 100d0
if nplanets eq 0 then k.derive=0

logk = parameter
logk.value = 10d0
logk.unit = 'm/s'
logk.description = 'RV semi-amplitude'
logk.latex = 'LOGK'
logk.label = 'logk'
logk.cgs = 100d0
logk.derive=0

period = parameter
period.unit = 'days'
period.description = 'Period'
period.latex = 'P'
period.label = 'Period'
if nplanets eq 0 then period.derive=0

a = parameter
a.unit = 'AU'
a.description = 'Semi-major axis'
a.latex = 'a'
a.label = 'a'
if nplanets eq 0 then a.derive=0

arsun = parameter
arsun.unit = '\rsun'
arsun.description = 'Semi-major axis in solar radii'
arsun.latex = 'arsun'
arsun.label = 'arsun'
arsun.cgs = 6.955d10
arsun.derive = 0

f0 = parameter
f0.value = 1d0
f0.description = 'Baseline flux'
f0.latex = 'F_0'
f0.label = 'f0'
f0.scale = 1d-2
if nplanets eq 0 then f0.derive = 0 $
else f0.fit = 1

detrend = parameter
detrend.value = 0d0
detrend.description = 'Detrending variable'
detrend.latex = ''
detrend.label = 'detrend'
detrend.scale = 0.1
detrend.fit = 1
detrend.derive = 1

beam = parameter
beam.value = 0d0
beam.description = 'Doppler Beaming amplitude'
beam.latex = 'A_B'
beam.unit = 'ppm'
beam.label = 'beam'
beam.scale = 100d0
beam.fit = 0
beam.derive = 0

bs = parameter
bs.description = 'Eclipse impact parameter'
bs.latex = 'b_S'
bs.label = 'bs'
if nplanets eq 0 then bs.derive = 0

ps = parameter
ps.description = 'A priori non-grazing eclipse prob'
ps.latex = 'P_S'
ps.label = 'ps'
if nplanets eq 0 then ps.derive = 0

psg = parameter
psg.description = 'A priori eclipse prob'
psg.latex = 'P_{S,G}'
psg.label = 'psg'
if nplanets eq 0 then psg.derive = 0

pt = parameter
pt.description = 'A priori non-grazing transit prob'
pt.latex = 'P_T'
pt.label = 'pt'
if nplanets eq 0 then pt.derive = 0

ptg = parameter
ptg.description = 'A priori transit prob'
ptg.latex = 'P_{T,G}'
ptg.label = 'ptg'
if nplanets eq 0 then ptg.derive = 0

mp = parameter
mp.unit = '\mj'
mp.description = 'Mass'
mp.latex = 'M_P'
mp.label = 'mp'
mp.cgs = 1.89813d30
if nplanets eq 0 then mp.derive = 0

mpearth = parameter
mpearth.unit = '\me'
mpearth.description = 'Mass'
mpearth.latex = 'M_P'
mpearth.label = 'mpearth'
mpearth.cgs = 5.97219d27
if nplanets eq 0 then mpearth.derive = 0

mpsun = parameter
mpsun.unit = '\msun'
mpsun.description = 'Mass'
mpsun.latex = 'M_P'
mpsun.label = 'mpsun'
mpsun.cgs = 1.9891d33
mpsun.derive = 0
mpsun.value = 0.001d0
if nplanets ne 0 then mpsun.fit = 1

msini = parameter
msini.unit = '\mj'
msini.description = 'Minimum mass'
msini.latex = 'M_P\sin i'
msini.label = 'msini'
msini.cgs = 1.89813d30
if nplanets eq 0 then msini.derive = 0

msiniearth = parameter
msiniearth.unit = '\me'
msiniearth.description = 'Minimum mass'
msiniearth.latex = 'M_P\sin i'
msiniearth.label = 'msiniearth'
msiniearth.cgs = 5.97219d27
if nplanets eq 0 then msiniearth.derive = 0

b = parameter
b.description = 'Transit impact parameter'
b.latex = 'b'
b.label = 'b'
if nplanets eq 0 then b.derive = 0

q = parameter
q.description = 'Mass ratio'
q.latex = 'M_P/M_*'
q.label = 'q'
if nplanets eq 0 then q.derive = 0

rp = parameter
rp.unit = '\rj'
rp.description = 'Radius'
rp.latex = 'R_P'
rp.label = 'rp'
rp.cgs = 7.1492d9
if nplanets eq 0 then rp.derive = 0

rpsun = parameter
rpsun.unit = '\rsun'
rpsun.description = 'Radius'
rpsun.latex = 'R_P'
rpsun.label = 'rpsun'
rpsun.cgs = 6.955d10
rpsun.derive = 0

rpearth = parameter
rpearth.unit = '\re'
rpearth.description = 'Radius'
rpearth.latex = 'R_P'
rpearth.label = 'rpearth'
rpearth.cgs = 6.3781d8
if nplanets eq 0 then rpearth.derive = 0

if keyword_set(earth) then begin
   ;; don't display jupiter unit parameters
   rp.derive = 0
   mp.derive = 0
   msini.derive = 0
endif else begin
   ;; don't display earth unit parameters
   rpearth.derive = 0
   mpearth.derive = 0
   msiniearth.derive = 0
endelse

rhop = parameter
rhop.unit = 'cgs'
rhop.description = 'Density'
rhop.latex = '\rho_P'
rhop.label = 'rhop'
if nplanets eq 0 then rhop.derive = 0

loggp = parameter
loggp.unit = 'cgs'
loggp.description = 'Surface gravity'
loggp.latex = 'logg_P'
loggp.label = 'loggp'
if nplanets eq 0 then loggp.derive = 0

teq = parameter
teq.unit = 'K'
teq.description = 'Equilibrium temp'
teq.latex = 'T_{\rm eq}'
teq.label = 'teq'
if nplanets eq 0 then teq.derive = 0

tcirc = parameter
tcirc.unit = 'Gyr'
tcirc.description = 'Tidal circ timescale'
tcirc.latex = '\tau_{\rm circ}'
tcirc.label = 'tcirc'
if nplanets eq 0 then tcirc.derive = 0

safronov = parameter
safronov.description = 'Safronov Number'
safronov.latex = '\Theta'
safronov.label = 'safronov'
if nplanets eq 0 then safronov.derive = 0

fave = parameter
fave.unit = '\fluxcgs'
fave.description = 'Incident Flux'
fave.latex = '\fave'
fave.label = 'fave'
fave.cgs = 1d-9
if nplanets eq 0 then fave.derive = 0

u1 = parameter
u1.description = 'Linear limb-darkening coeff'
u1.latex = 'u_{1}'
u1.label = 'u1'
u1.scale = 0.15d0
if nplanets eq 0 then u1.derive = 0 $
else u1.fit = 1

u2 = parameter
u2.description = 'Quadratic limb-darkening coeff'
u2.latex = 'u_{2}'
u2.label = 'u2'
u2.scale = 0.15d0
if nplanets eq 0 then u2.derive = 0 $
else u2.fit = 1

u3 = parameter
u3.description = 'Non-linear limb-darkening coeff'
u3.latex = 'u_{3}'
u3.label = 'u3'
u3.scale = 0.15d0
u3.derive = 0

u4 = parameter
u4.description = 'non-linear limb-darkening coeff'
u4.latex = 'u_{4}'
u4.label = 'u4'
u4.scale = 0.15d0
u4.derive = 0

thermal = parameter
thermal.description = 'Thermal emission from the planet'
thermal.latex = 'A_T'
thermal.label = 'thermal'
thermal.scale = 1d4
thermal.unit = 'ppm'
thermal.fit = 0
thermal.derive = 0

reflect = parameter
reflect.description = 'Reflection from the planet'
reflect.latex = 'A_R'
reflect.label = 'reflect'
reflect.scale = 1d4
reflect.unit = 'ppm'
reflect.fit = 0
reflect.derive = 0

phaseshift = parameter
phaseshift.description = 'Phase shift of reflection curve'
phaseshift.latex = '\phi'
phaseshift.label = 'phaseshift'
phaseshift.scale = 90d0
phaseshift.unit = 'deg'
phaseshift.fit = 0
phaseshift.derive = 0

ellipsoidal = parameter
ellipsoidal.description = 'Ellipsoidal variation of the star'
ellipsoidal.latex = 'A_E'
ellipsoidal.label = 'ellipsoidal'
ellipsoidal.scale = 1d4
ellipsoidal.unit = 'ppm'
ellipsoidal.fit = 0
ellipsoidal.derive = 0

dilute = parameter
dilute.description = 'Dilution from neighboring stars'
dilute.latex = 'A_D'
dilute.label = 'dilute'
dilute.scale = 1d-2
dilute.fit = 0
dilute.derive = 0

mag = parameter
mag.description = 'Apparent Magnitude'
mag.latex = 'mag'
mag.label = 'mag'
mag.scale = 0.1
mag.fit = 0
mag.derive = 0

ttv = parameter
ttv.description = 'Transit Timing Variation'
ttv.latex = 'TTV'
ttv.label = 'ttv'
ttv.unit = 'days'
ttv.scale = 0.02 ;; ~30 minutes
junk = where(fittran, nfittran)
ttv.derive=0

tiv = parameter
tiv.description = 'Transit Inclination Variation'
tiv.latex = 'TiV'
tiv.label = 'tiv'
tiv.unit = 'Radians'
tiv.scale = 0.01*!dpi/180d0
tiv.derive=0

tdeltav = parameter
tdeltav.description = 'Transit Depth Variation'
tdeltav.latex = 'T\delta V'
tdeltav.label = 'tdeltav'
tdeltav.scale = 0.01
tdeltav.derive=0

jitter = parameter
jitter.description = 'RV Jitter'
jitter.latex = '\sigma_J'
jitter.label = 'jitter'
jitter.unit = 'm/s'
jitter.value = 0d0
jitter.scale = 1d0
jitter.fit = 0
if nplanets eq 0 then jitter.derive = 0

jittervar = parameter
jittervar.description = 'RV Jitter Variance'
jittervar.latex = '\sigma_J^2'
jittervar.label = 'jittervar'
jittervar.value = 0d0
jittervar.scale = 1d0
if nplanets eq 0 then jitter.derive = 0 $
else jittervar.fit=1

variance = parameter
variance.description = 'Added Variance'
variance.latex = '\sigma^{2}'
variance.label = 'variance'
variance.value = 0d0
variance.scale = 1d0
if nplanets eq 0 then variance.derive = 0 $
else variance.fit=1

sigma_r = parameter
sigma_r.description = 'Red noise'
sigma_r.latex = '\sigma_r'
sigma_r.label = 'sigma_r'
sigma_r.value = 0d0
sigma_r.scale = 1d0
sigma_r.derive = 0
sigma_r.fit=0

errscale = parameter
errscale.description = 'SED photometry error scaling'
errscale.latex = '\sigma_{SED}'
errscale.label = 'errscale'
errscale.value = 1d0
errscale.scale = 10d0
errscale.fit = 0
errscale.derive = 0

sperrscale = parameter
sperrscale.description = 'SED spectrophotometry error scaling'
sperrscale.latex = '\sigma_{SP}'
sperrscale.label = 'sperrscale'
sperrscale.value = 1d0
sperrscale.scale = 10d0
sperrscale.fit = 0
sperrscale.derive = 0

;; Create the structures -- The order here dictates the order in the
;;                          output table.

;; for each star
columnlabels = ''
star = create_struct(mstar.label,mstar,$
                     rstar.label,rstar,$
                     rstarsed.label,rstarsed,$
                     lstar.label,lstar,$
                     fbol.label,fbol,$
                     rhostar.label,rhostar,$
                     logg.label,logg,$
                     teff.label,teff,$
                     teffsed.label,teffsed,$
                     feh.label,feh,$
                     fehsed.label,fehsed,$
                     initfeh.label,initfeh,$
                     age.label,age,$
                     eep.label,eep,$
                     logmstar.label,logmstar,$
                     absks.label,absks,$
                     appks.label,appks,$
                     vsini.label,vsini,$
                     vline.label,vline,$                    
                     Av.label,Av,$
                     alpha.label,alpha,$
;                     Ma.label,Ma,$
;                     Mv.label,Mv,$
                     errscale.label,errscale,$
                     ra.label,ra,$       ;; astrometry
                     dec.label,dec,$     ;; astrometry
                     pmra.label,pmra,$   ;; astrometry
                     pmdec.label,pmdec,$ ;; astrometry
                     rvabs.label,rvabs,$ ;; astrometry
                     parallax.label,parallax,$
                     distance.label,distance,$
                     slope.label,slope,$
                     quad.label,quad,$
                     'rootlabel','Stellar Parameters:',$
                     'columnlabels',columnlabels,$
                     'label','')
       
;; Unit constants (converted to cgs)
;; As defined by IAU resolutions B2, B3
;; https://www.iau.org/static/resolutions/IAU2012_English.pdf
;; https://arxiv.org/abs/1510.07674, 
;; https://arxiv.org/abs/1507.07956, Table 1
constants = mkconstants()

ndata = 0L


;if n_elements(nvalues) ne 0 then stop

planet = create_struct($
         period.label,period,$    ;; fundamental (most interesting) parameters
         rp.label,rp,$
         rpearth.label,rpearth,$
         mp.label,mp,$
         mpsun.label,mpsun,$
         logmp.label,logmp,$
         mpearth.label,mpearth,$
         tc.label,tc,$
         tt.label,tt,$
         t0.label,t0,$
         a.label,a,$              
         inc.label,inc,$
         ideg.label,ideg,$
         e.label,e,$
         omega.label,omega,$
         omegadeg.label,omegadeg,$
         lsinw.label,lsinw,$
         lsinw2.label,lsinw2,$
         lcosw.label,lcosw,$
         bigomega.label,bigomega,$ ;; for astrometry
         bigomegadeg.label,bigomegadeg,$
         lsinbigomega.label,lsinbigomega,$
         lcosbigomega.label,lcosbigomega,$
         teq.label,teq,$
         tcirc.label, tcirc,$
         K.label,k,$              ;; RV parameters
         logK.label,logk,$              ;; RV parameters
         p.label,p,$              ;; Primary Transit parameters
         ar.label,ar,$
         delta.label,delta)

;; compute a depth for each observed band
for i=0L, nband-1 do begin
   depth = parameter
   prettyname = prettybands[(where(bands[i] eq allowedbands))[0]]
   depth.description = 'Transit depth in ' + prettyname
   depth.latex = '\delta_{\rm ' + prettyname + '}'
   depth.label = 'depth_' + bands[i]
   depth.unit = 'frac'
   if nplanets eq 0 then depth.derive = 0
   planet = create_struct(planet, depth.label, depth)
endfor

planet = create_struct($
         planet,$
         tau.label,tau,$
         t14.label,t14,$
         tfwhm.label,tfwhm,$
         b.label,b,$
         cosi.label,cosi,$
         bs.label,bs,$          ;; secondary eclipse parameters
         taus.label,taus,$
         t14s.label,t14s,$
         tfwhms.label,tfwhms,$
;         eclipsedepth36.label,eclipsedepth36,$
;         eclipsedepth45.label,eclipsedepth45,$
         eclipsedepth25.label,eclipsedepth25,$
         eclipsedepth50.label,eclipsedepth50,$
         eclipsedepth75.label,eclipsedepth75,$
         rhop.label,rhop,$      ;; less useful parameters            
         rpsun.label,rpsun,$
         logP.label,logp,$  
         loggp.label,loggp,$
         lambda.label,lambda,$
         lambdadeg.label,lambdadeg,$
         lsinlambda.label,lsinlambda,$
         lcoslambda.label,lcoslambda,$
         safronov.label,safronov,$
         fave.label,fave,$
         tp.label,tp,$
         ts.label,ts,$
         ta.label,ta,$
         td.label,td,$
         phase.label,phase,$
         vcve.label,vcve,$
         chord.label,chord,$
         sign.label,sign,$
         sign2.label,sign2,$
         ecosw.label,ecosw,$
         esinw.label,esinw,$
         secosw.label,secosw,$
         sesinw.label,sesinw,$
         qecosw.label,qecosw,$
         qesinw.label,qesinw,$
         msini.label,msini,$
         msiniearth.label,msiniearth,$
         q.label,q,$
         arsun.label,arsun,$
         dr.label,dr,$
         pt.label,pt,$
         ptg.label,ptg,$
         ps.label,ps,$                 
         psg.label,psg,$     
         beam.label,beam,$     ;; other
         'starndx',0L,$
         'fittran',fittran[0],$        ;; booleans
         'fitrv',fitrv[0],$
         'chen',chen[0],$
         'i180',i180[0],$
         'rossiter',rossiter[0],$
         'fitdt',fitdt[0],$
         'rootlabel','Planetary Parameters:',$
         'label','')

specphot = create_struct(sperrscale.label, sperrscale) ;; spectrophotometry error scaling

;; for each wavelength
band = create_struct(u1.label,u1,$ ;; linear limb darkening
                     u2.label,u2,$ ;; quadratic limb darkening
                     u3.label,u3,$ ;; 1st non-linear limb darkening
                     u4.label,u4,$ ;; 2nd non-linear limb darkening
                     thermal.label,thermal,$ ;; thermal emission
                     ;dilute.label,dilute,$   ;; dilution
                     reflect.label,reflect,$ ;; reflection
                     phaseshift.label,phaseshift,$ ;; reflection
                     ellipsoidal.label,ellipsoidal,$
                     eclipsedepth.label,eclipsedepth,$
                     mag.label,mag,$
                     phottobary.label,phottobary,$
                     'starndx',0L,$
                     'name','',$
                     'rootlabel','Wavelength Parameters:',$
                     'label','')

;; for each telescope
telescope = create_struct(gamma.label,gamma,$
                          jitter.label,jitter,$
                          jittervar.label,jittervar,$
                          'rvptrs', ptr_new(),$
                          'detrend',ptr_new(/allocate_heap),$ ;; array of detrending parameters
                          'name','',$
                          'chi2',0L,$
                          'rootlabel','Telescope Parameters:',$
                          'label','')

if ntel le 0 then begin
   telescope.jittervar.fit = 0
   telescope.jittervar.derive = 0
   telescope.jitter.derive = 0
endif

;; for each transit
transit = create_struct(variance.label,variance,$ ;; jitter
                        sigma_r.label,sigma_r,$ ;; Red noise
                        ttv.label,ttv,$           ;; Transit Timing Variation
                        tiv.label,tiv,$ ;; Transit inclination variation
                        tdeltav.label,tdeltav,$ ;; Transit depth variation
                        dilute.label, dilute, $ ;; Transit depth variation
                        f0.label,f0,$ ;; normalization
                        'claret', 0B,$
                        'transitptrs',ptr_new(),$ ;; Data
                        'detrend',ptr_new(/allocate_heap),$ ;; array of detrending parameters
                        'bandndx',0L,$
                        'exptime',0d0,$
                        'ninterp',1d0,$
                        'name','',$
                        'epoch',dblarr(nplanets>1) + !values.d_nan,$
                        'pndx',0L,$ ;; index to which planet this corresponds to (-1=>all)
                        'bitmask',0ULL,$ ;; a bitmask for which planet(s) this transit corresponds to
                        'chi2',0L,$
                        'rootlabel','Transit Parameters:',$
                        'rejectflatmodel',0B,$
                        'fitspline',0B,$
                        'splinespace',0.75d0,$
                        'label','') 

doptom = create_struct('dtptrs',ptr_new(),$
                       'rootlabel','Doppler Tomography Parameters:',$
                       'label','',$
                       dtscale.label,dtscale) ;,$
;                       'lambdarange',0,$
;                       'tel','',$
;                       'night','',$
;                       'planetndx')

astrom = create_struct('astromptrs',ptr_new(),$
                       astromscale.label,astromscale,$
                       raoffset.label,raoffset,$
                       decoffset.label,decoffset,$
                       'starndx',0L,$
                       'bandndx',0L,$
                       'rootlabel','Astrometry Parameters:',$
                       'label','',$
                       'epoch',0L)

;; a stellar system has a star, planets, observed bands, observed
;; transits, priors, and global options
ss = create_struct('star',replicate(star,nstars>1),$
                   'planet',replicate(planet,nplanets > 1),$
                   'band',replicate(band,nband > 1),$
                   'telescope',replicate(telescope,ntel > 1),$
                   'transit',replicate(transit,ntran > 1),$
                   'doptom',replicate(doptom,ndt>1),$
                   'astrom',replicate(astrom,nastrom>1),$
                   'specphot',replicate(specphot,nspecfiles>1),$
                   'epochs',dblarr(ntran>1,nplanets>1),$
                   'constants',constants,$
                   'tofit',ptr_new(1),$
                   'priors',ptr_new(1),$
                   'debug',keyword_set(debug),$
                   'verbose',keyword_set(verbose),$
                   'tides',keyword_set(tides),$
                   'ntel',ntel,$
                   'rvepoch',0d0,$
                   'ntran',ntran,$
                   'nastrom',nastrom,$
                   'nband',nband,$
                   'nstars',nstars,$
                   'ndt',ndt,$
                   'nplanets',nplanets,$
                   'ndata',ndata,$
                   'mist', mist,$
                   'parsec', parsec,$
                   'yy', yy,$
                   'torres', torres,$
                   'mannrad', mannrad,$
                   'mannmass', mannmass,$
                   'fbolsedfloor', fbolsedfloor,$
                   'teffsedfloor', teffsedfloor,$
                   'fehsedfloor', fehsedfloor,$
                   'teffemfloor', teffemfloor,$
                   'fehemfloor', fehemfloor,$
                   'rstaremfloor', rstaremfloor,$
                   'ageemfloor', ageemfloor,$
                   'avprior', avprior, $
                   'diluted',diluted,$
                   'dilutebandndx',ptr_new(1),$;[-1L,1,1],$
                   'dilutestarndx',dilutestarndx,$
                   'oned', keyword_set(oned),$
                   'ttvs', ttvs,$
                   'tivs', tivs,$
                   'tdeltavs', tdeltavs,$
                   'alloworbitcrossing', keyword_set(alloworbitcrossing),$
                   'nsteps',nsteps,$                   
                   'npars',0L,$
                   'burnndx',0L,$
                   'nchains',1L,$
                   'goodchains',ptr_new(1),$
                   'amoeba',0L,$
                   'logname','',$
                   'chi2',ptr_new(1),$
                   'fluxfile',fluxfile,$
                   'mistsedfile',mistsedfile,$
                   'sedfile',sedfile,$
                   'specphotpath',specphotpath,$
                   'nspecfiles',nspecfiles,$
                   ;; metadata to be able to restart fit
                   'circular', circular,$
                   'fitrv',fitrv,$
                   'fittran',fittran,$
                   'fitdt',fitdt,$
                   'rossiter',rossiter,$
                   'fitlogmp',fitlogmp,$
                   'rejectflatmodel',rejectflatmodel,$
                   'longcadence',longcadence,$
                   'tranpath',tranpath,$
                   'rvpath',rvpath,$
                   'astrompath',astrompath,$
                   'fitslope',keyword_set(fitslope),$
                   'fitquad',keyword_set(fitquad),$
                   'fitspline',fitspline,$
                   'splinespace',splinespace,$
                   'fitwavelet',fitwavelet,$
                   'ninterp',ninterp,$
                   'exptime',exptime,$
                   'novcve',novcve,$
                   'nochord',nochord,$
                   'fitsign',fitsign,$
                   'chi2func',chi2func,$
                   'fittt',fittt,$
                   'i180',i180,$
                   'chen',chen,$
                   'fitdilute',fitdilute,$
                   'fitthermal',fitthermal,$
                   'fitreflect',fitreflect,$
                   'fitphase',fitphase,$
                   'fitellip',fitellip,$
                   'fitbeam',fitbeam,$
                   'derivebeam',derivebeam,$
                   'planetorder',lindgen(nplanets > 1),$
                   'dtpath',dtpath,$
                   'delay',delay,$
                   'prefix',prefix, $
                   'earth',keyword_set(earth),$
                   'priorfile',priorfile,$
                   'transitrange',transitrange,$
                   'rvrange',rvrange,$
                   'sedrange',sedrange,$
                   'emrange',emrange)

;)

;; three different ways to do the SED model
if file_test(mistsedfile) or file_test(sedfile) or file_test(fluxfile) then begin

   ;; overwrite the common block, in case it's been called
   ;; before then updated (without exiting IDL)
   ;; the chi2 doesn't matter here, use solar values
   
   if file_test(mistsedfile) then begin
;      common BC_block, bcarrays, teffgrid, logggrid, fehgrid, avgrid, sedbands, mags, errs, filterprops, blend
      sedchi2 = mistmultised(replicate(6000d0,nstars), replicate(4.41d0,nstars), replicate(0d0,nstars), $
                             replicate(0d0,nstars), replicate(10d0,nstars), replicate(1d0,nstars), $
                             replicate(1d0,nstars), mistsedfile, /redo,blend0=blend)
      ss.mistsedfile = mistsedfile
      ss.ndata += n_elements(mags) +2d0 ;; two more because of links between rstar and rstarsed, teff and teffsed
   endif else if file_test(sedfile) then begin
      sedchi2 = exofast_multised(replicate(6000d0,nstars), replicate(4.41d0,nstars), replicate(0d0,nstars), $
                                 replicate(0d0,nstars), replicate(10d0,nstars), replicate(1d0,nstars), $
                                 replicate(1d0,nstars), sedfile, /redo, specphotpath=specphotpath,$
                                 blend0=blend,rstar=replicate(1d0,nstars),sperrscale=ss.specphot.sperrscale.value)
      ss.sedfile = sedfile
      readcol, sedfile, junk, format='a', comment='#', /silent
      ndata += n_elements(junk) + 2 ;; two more because of links between rstar and rstarsed, teff and teffsed
      for i=0L, ss.nspecfiles-1 do begin
         ss.specphot.sperrscale.fit = 1
         ss.specphot.sperrscale.derive = 1
      endfor
   endif else begin
      printandlog, 'WARNING: FLUXFILE has been deprecated. MISTSEDFILE should be used instead.', logname
      printandlog, 'NOTE: When using MISTSEDFILE, the C3K atmosphere is not computed directly. The plotted NextGen atmosphere is for aesthetics only.', logname

      sedchi2 = exofast_sed(fluxfile, 6000d0,1d0,0d0,10d0,logg=4.41d0,met=0d0,alpha=0d0,/redo)
      ss.fluxfile = fluxfile
      readcol, fluxfile, junk, format='a', comment='#', /silent
      ndata += n_elements(junk) + 2 ;; two more because of links between rstar and rstarsed, teff and teffsed
   endelse
   
   starsused = where(total(blend,1) ne 0,nused)
   for i=0L, nused-1 do begin
      
      ss.star[starsused[i]].errscale.fit = 1
      ss.star[starsused[i]].errscale.derive = 1
      ss.star[starsused[i]].distance.derive = 1
      ss.star[starsused[i]].distance.fit = 1
      ss.star[starsused[i]].fbol.derive = 1
      ss.star[starsused[i]].parallax.derive = 1
      ss.star[starsused[i]].av.fit = 1
      ss.star[starsused[i]].av.derive = 1
      
      if teffsedfloor ne 0d0 then begin
         ss.star[starsused[i]].teffsed.fit = 1
         ss.star[starsused[i]].teffsed.derive = 1
      endif
      if fbolsedfloor ne 0d0 then begin
         ss.star[starsused[i]].rstarsed.fit = 1
         ss.star[starsused[i]].rstarsed.derive = 1
      endif
      if fehsedfloor ne 0d0 then begin
         ss.star[starsused[i]].fehsed.fit = 1
         ss.star[starsused[i]].fehsed.derive = 1
      endif
   endfor
endif else begin
   if mistsedfile ne '' then begin
      printandlog, 'Could not find ' + mistsedfile, logname
      return, -1
   endif else if sedfile ne '' then begin
      printandlog, 'Could not find ' + sedfile, logname
      return, -1
   endif else if fluxfile ne '' then begin
      printandlog, 'Could not find ' + fluxfile, logname
      return, -1
   endif
endelse

if n_elements(logname) eq 1 then ss.logname=logname

for i=0L, nstars-1 do begin
   if ss.mist[i] or ss.parsec[i] then begin
      ss.star[i].initfeh.fit=1
      ss.star[i].initfeh.derive=1
      ss.star[i].age.derive=1
      ss.star[i].age.fit=1
      ss.star[i].eep.fit=1
      ss.star[i].eep.derive=1
   endif else begin
      ss.star[i].initfeh.fit=0
      ss.star[i].initfeh.derive=0
      ss.star[i].eep.fit=0
      ss.star[i].eep.derive=0      
   endelse

   if ss.mannrad[i] or ss.mannmass[i] then begin
      ss.star[i].distance.fit = 1
      ss.star[i].distance.derive = 1
      ss.star[i].parallax.derive = 1
      ss.star[i].appks.fit = 1
      ss.star[i].appks.derive = 1
      ss.star[i].absks.derive = 1
   endif

   if ss.yy[i] then begin
      ss.star[i].age.fit=1
      ss.star[i].age.derive=1
   endif
   ss.star[i].label = strtrim(i,2)

   ;; these are global parameters, not one per star
   ;; Should move to SS structure, but may require crawling routine modifications?
   if i gt 0 then begin
      ss.star[i].slope.fit = 0
      ss.star[i].slope.derive = 0
      ss.star[i].quad.fit = 0
      ss.star[i].quad.derive = 0
      ss.star[i].errscale.fit = 0
      ss.star[i].errscale.derive = 0
   endif

endfor

if ndt gt 0 then begin
   ss.doptom[*].dtptrs = ptrarr(ndt,/allocate_heap)
   if ~keyword_set(silent) then printandlog, 'The index for each DT file is',logname
   for i=0, ndt-1 do begin
      *(ss.doptom[i].dtptrs) = exofast_readdt(dtfiles[i])
      if ~keyword_set(silent) then printandlog, string(i,dtfiles[i],format='(i2,x,a)'),logname
   endfor
   if ~keyword_set(silent) then printandlog, '', logname
endif

;; for each planet 
;; populate the planet fitting parameters
;; planetary labels, a bit optimistic...
plabels = ['b','c','d','e','f','g','h','i','j','k','l','m','n',$
           'o','p','q','r','s','t','u','v','w','x','y','z'] 

for i=0, nplanets-1 do begin
   ss.planet[i].label = plabels[i]
   ss.planet[i].starndx = starndx[i]

   ;; circular orbit, don't fit e or omega
   if circular[i] then begin
      ss.planet[i].qesinw.fit = 0
      ss.planet[i].qecosw.fit = 0
      ss.planet[i].sesinw.fit = 0
      ss.planet[i].secosw.fit = 0
      ss.planet[i].e.derive=0
      ss.planet[i].omegadeg.derive=0
      ss.planet[i].esinw.derive=0
      ss.planet[i].ecosw.derive=0
   
      ;; don't derive secondary parameters (same as primary if circular)
      ss.planet[i].taus.derive = 0
      ss.planet[i].t14s.derive = 0
      ss.planet[i].tfwhms.derive = 0
      ss.planet[i].bs.derive = 0
      ss.planet[i].ps.derive = 0
      ss.planet[i].psg.derive = 0

   endif
   
   ;; if rossiter is done, fit lambda (for each planet) and vsini (for the star)
   ;; if DT is done, fit lambda (for each planet) and vsini and vline (for the star)
   if rossiter[i] or fitdt[i] then begin
      ss.star[ss.planet[i].starndx].vsini.fit = 1
      ss.star[ss.planet[i].starndx].vsini.derive = 1
      ss.planet[i].lsinlambda.fit = 1
      ss.planet[i].lcoslambda.fit = 1
      ss.planet[i].lambdadeg.derive = 1
      if fitdt[i] then begin
         ss.planet[i].fitdt = 1B
         ss.star[ss.planet[i].starndx].vline.fit = 1
         ss.star[ss.planet[i].starndx].vline.derive = 1
         ss.doptom[*].dtscale.fit = 1
         ss.doptom[*].dtscale.derive = 1
      endif
      if rossiter[i] then ss.planet[i].rossiter = 1B
   endif

   ss.planet[i].fittran = fittran[i]
   ss.planet[i].fitrv = fitrv[i]
   ss.planet[i].chen = chen[i]
   if ss.planet[i].fittran then ss.planet[i].tt.derive = 1B $
   else ss.planet[i].tt.derive = 0B

   if ss.planet[i].fittran and fittt[i] then begin
      ss.planet[i].tt.fit = 1B
      ss.planet[i].tc.fit = 0B
   endif

   ;; fit in Vc/Ve, lcosw, lsinw, and sign if only fitting transit (far more efficient)
   ;; fit in tfwhm0, lcosw, lsinw, and sign if only fitting transit (far more efficient)
   if ~ss.planet[i].fitrv and ss.planet[i].fittran and ~circular[i] and ~novcve[i] then begin
      ss.planet[i].secosw.fit = 0
      ss.planet[i].sesinw.fit = 0
      ss.planet[i].lsinw.fit = 0
      ss.planet[i].lsinw2.fit = 0
      ss.planet[i].lcosw.fit = 0

      ;; reparameterize secosw, sesinw => vcve, lsinw, lcosw
      ss.planet[i].lsinw.fit = 1
      ss.planet[i].lsinw2.fit = 0
      ss.planet[i].lcosw.fit = 1
      ss.planet[i].vcve.fit = 1
      ss.planet[i].vcve.derive = 1
      
      ;; use sign to choose between + and - solutions to the quadratic
      ;; to solve for e (otherwise, use lsinw^2 + lcosw^2 = L)
      if fitsign[i] then begin
         ss.planet[i].sign.fit = 1
         ss.planet[i].sign.derive = 1
;         ss.planet[i].sign2.fit = 1
;         ss.planet[i].sign2.derive = 1
      endif

      ;; reparameterize cosi as chord
      if ~nochord[i] then begin
         ss.planet[i].cosi.fit = 0
         ss.planet[i].chord.fit = 1         
         ss.planet[i].chord.derive = 1
      endif
   endif

   ;; we can marginalize over these parameters 
   ;; even if a transit is not fit.
   ;; And with the Chen & Kipping relation, we can get a decent value
   ;; on the radius and density
   if not fittran[i] and chen[i] then ss.planet[i].cosi.scale = 1d0

   ;; Chen & Kipping prior can be used to fit planetary radius 
   ;; and we can marginalize over cosi
   if not chen[i] then begin
      if not fittran[i] and nastrom eq 0 then begin
         ss.planet[i].cosi.fit = 0
         ss.planet[i].p.fit = 0
         ss.planet[i].cosi.derive = 0
        
         ss.planet[i].p.derive = 0
         ss.planet[i].ideg.derive = 0
         ss.planet[i].delta.derive = 0
         ;;ss.planet[i].depth.derive = 0
         ss.planet[i].mp.derive = 0
         ss.planet[i].mpearth.derive = 0
         
         ;; transit derived pars
         ss.planet[i].b.derive = 0
         ss.planet[i].tfwhm.derive = 0
         ss.planet[i].t14.derive = 0
         ss.planet[i].tau.derive = 0
         ss.planet[i].bs.derive = 0
         
         ;; eclipse derived pars
         ss.planet[i].tfwhms.derive = 0
         ss.planet[i].t14s.derive = 0
         ss.planet[i].taus.derive = 0
         
         ;; requires mass and radius
         ss.planet[i].rhop.derive = 0
         ss.planet[i].loggp.derive = 0
         ss.planet[i].safronov.derive = 0
      endif

   endif

   if fitbeam[i] then begin
      ss.planet[i].beam.fit = 1B
      ss.planet[i].beam.derive = 1B
   endif else if derivebeam[i] then ss.planet[i].beam.derive = 1B

   if i180[i] or ss.nastrom gt 0 then ss.planet[i].i180 = 1

   if fitlogmp[i] then begin
      ss.planet[i].logmp.fit = 1
      ss.planet[i].mpsun.fit = 0
   endif

endfor

for i=0, nband-1 do begin
   ss.band[i].name = bands[i]
   ss.band[i].label = prettybands[(where(bands[i] eq allowedbands))[0]]

   ldcoeffs = quadld(ss.star[0].logg.value, ss.star[0].teff.value, ss.star[0].feh.value, bands[i])
   if finite(ldcoeffs[0]) then ss.band[i].u1.value = ldcoeffs[0] $
   else ss.band[i].u1.value = 0d0
   if finite(ldcoeffs[1]) then ss.band[i].u2.value = ldcoeffs[1] $
   else ss.band[i].u2.value = 0d0

   match = where(fitthermal eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].thermal.fit = 1B
      ss.band[i].thermal.derive = 1B
      ss.band[i].eclipsedepth.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting thermal emission for " + ss.band[i].name + " band",logname
   endif

   match = where(fitreflect eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].reflect.fit = 1B
      ss.band[i].reflect.derive = 1B
      ss.band[i].eclipsedepth.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting reflected light for " + ss.band[i].name + " band", logname
   endif

;   match = where(fitphase eq ss.band[i].name)
;   if match[0] ne -1 then begin
;      ss.band[i].phaseshift.fit = 1B
;      ss.band[i].phaseshift.derive = 1B
;      if ~keyword_set(silent) then printandlog, "Fitting phase offset for " + ss.band[i].name + " band", logname
;   endif

;   match = where(fitdilute eq ss.band[i].name)
;   if match[0] ne -1 then begin
;      ss.band[i].dilute.fit = 1B
;      ss.band[i].dilute.derive = 1B
;      if ~keyword_set(silent) then printandlog, "Fitting dilution for " + ss.band[i].name + " band", logname
;   endif

   match = where(fitellip eq ss.band[i].name)
   if match[0] ne -1 then begin
      ss.band[i].ellipsoidal.fit = 1B
      ss.band[i].ellipsoidal.derive = 1B
      if ~keyword_set(silent) then printandlog, "Fitting ellipsoidal for " + ss.band[i].name + " band", logname
   endif

endfor

if nband eq 0 then begin
   ss.band[0].u1.fit=0B
   ss.band[0].u2.fit=0B
   ss.band[0].u1.derive=0B
   ss.band[0].u2.derive=0B
endif

minallbjd = !values.d_infinity
maxallbjd = -!values.d_infinity
dilutebandndx = [-1]

;; read in the transit files
if ntran gt 0 then begin

   if ~keyword_set(silent) then printandlog, 'The index for each transit is',logname
   ss.transit[*].transitptrs = ptrarr(ntran,/allocate_heap)
   for i=0, ntran-1 do begin
      if ~keyword_set(silent) then printandlog, string(i,tranfiles[i],format='(i2,x,a)'),logname
      *(ss.transit[i].transitptrs) = readtran(tranfiles[i], detrendpar=detrend, nplanets=nplanets,skipallowed=noclaret[i])

;      ;; create an array of detrending variables 
;      ;; (one for each extra column in the transit file)
      nadd = (*(ss.transit[i].transitptrs)).nadd
;      if nadd ge 1 then *(ss.transit[i].detrendadd) = replicate(detrendadd,nadd)
      nmult = (*(ss.transit[i].transitptrs)).nmult
;      if nmult ge 1 then *(ss.transit[i].detrendmult) = replicate(detrendmult,nmult)

      ss.ndata += n_elements((*(ss.transit[i].transitptrs)).bjd)*(1L+nadd+nmult)

      minbjd = min((*(ss.transit[i].transitptrs)).bjd,max=maxbjd)
      if minbjd lt minallbjd then minallbjd = minbjd
      if maxbjd gt maxallbjd then maxallbjd = maxbjd


      ss.transit[i].exptime = exptime[i]
      ss.transit[i].ninterp = ninterp[i]
      
      band = (*(ss.transit[i].transitptrs)).band
      ss.transit[i].bandndx = where(ss.band[*].name eq band)
      ss.transit[i].label = (*(ss.transit[i].transitptrs)).label
      ss.transit[i].rejectflatmodel = rejectflatmodel[i]

      if total(diluted[i,*]) then begin
         dilutebandndx = [dilutebandndx,ss.transit.bandndx]
         ss.transit[i].dilute.fit = 1B
         ss.transit[i].dilute.derive = 1B
      endif

      if fitdilute[i] then begin
         ss.transit[i].dilute.fit = 1B
         ss.transit[i].dilute.derive = 1B
      endif

      ss.transit[i].claret = ~keyword_set(noclaret[i])
      ss.transit[i].fitspline = fitspline[i]
      ss.transit[i].splinespace = splinespace[i]
;      ;; F0 and the spline are totally degenerate. 
;      ;; Don't fit (or derive) F0 if flattening with a spline
;      if ss.transit[i].fitspline then begin
;         ss.transit[i].f0.fit = 0B
;         ss.transit[i].f0.derive = 0B
;      endif         

      if fitwavelet[i] then begin
         ss.transit[i].variance.fit = 0B
         ss.transit[i].variance.derive = 0B
         ss.transit[i].sigma_r.fit = 1B
         ss.transit[i].sigma_r.derive = 1B
      endif

      for j=0L, ss.nplanets-1L do begin
         if keyword_set(ttvs[i,j]) then begin
            ss.transit[i].ttv.fit = 1B
            ss.transit[i].ttv.derive = 1B
         endif
         if keyword_set(tivs[i,j]) then begin
            ss.transit[i].tiv.fit = 1B
            ss.transit[i].tiv.derive = 1B
         endif
         if keyword_set(tdeltavs[i,j]) then begin
            ss.transit[i].tdeltav.fit = 1B
            ss.transit[i].tdeltav.derive = 1B
         endif
      endfor

   endfor

   if n_elements(dilutebandndx) gt 1 then begin
      dilutebandndx = dilutebandndx[1:*]
      ;; remove duplicates, assign to structure
      sorted = sort(dilutebandndx)
      unique = uniq(dilutebandndx[sorted])
      dilutebandndx = dilutebandndx[sorted[unique]]
   endif
   
   if ~keyword_set(silent) then printandlog, '',logname
endif else begin
   ss.transit[0].f0.fit = 0
   ss.transit[0].f0.derive = 0
   ss.transit[0].variance.fit = 0
   ss.transit[0].variance.derive = 0
endelse
*ss.dilutebandndx = dilutebandndx

;; read in the RV files
if ntel gt 0 then begin
   ss.telescope[*].rvptrs = ptrarr(ntel,/allocate_heap)
   if ~keyword_set(silent) then printandlog, "The index for each RV data set is",logname
   maxpoints = 0
   detrend.scale = 1d0
   for i=0, ntel-1 do begin
      if ~keyword_set(silent) then printandlog, string(i,rvfiles[i],format='(i2,x,a)'),logname
      *(ss.telescope[i].rvptrs) = readrv_detrend(rvfiles[i], detrendpar=detrend)
      ss.telescope[i].label = (*(ss.telescope[i].rvptrs)).label

      minbjd = min((*(ss.telescope[i].rvptrs)).bjd,max=maxbjd)
      if minbjd lt minallbjd then minallbjd = minbjd
      if maxbjd gt maxallbjd then maxallbjd = maxbjd

      ;; create an array of detrending variables 
      ;; (one for each extra column in the rv file)
      nadd = (*(ss.telescope[i].rvptrs)).nadd
      nmult = (*(ss.telescope[i].rvptrs)).nmult
      ss.ndata += n_elements((*(ss.telescope[i].rvptrs)).bjd)*(1L+nadd+nmult)

   endfor
   if ~keyword_set(silent) then printandlog, '', logname
endif else begin
   ss.telescope[*].rvptrs = ptr_new(/allocate_heap)
endelse

for i=0, ss.ntel-1 do begin
   rv = *(ss.telescope[i].rvptrs)
   ss.telescope[i].gamma.value = mean(rv.rv)

;   print, (*ss.telescope[i].rvptrs).planet, (*ss.telescope[i].rvptrs).planet eq -1
;stop

   if (*ss.telescope[i].rvptrs).planet eq -1 then begin
      if n_elements(alltime) eq 0 then begin
         alltime = rv.bjd
         allrv = rv.rv-ss.telescope[i].gamma.value
      endif else begin
         alltime = [alltime,rv.bjd]
         allrv = [allrv,rv.rv-ss.telescope[i].gamma.value]
      endelse
   endif
endfor

;; take a rough stab at gamma, slope, quad, and K
if ss.ntel gt 0 and n_elements(rvepoch) eq 0 then ss.rvepoch = (min(alltime) + max(alltime))/2d0
if ss.star[0].quad.fit then begin
   coeffs = poly_fit(alltime-ss.rvepoch, allrv, 2, yfit=yfit)
   allrv -= yfit
   ss.star[0].quad.value = coeffs[2]
   ss.star[0].slope.value = coeffs[1]
   ss.star[0].quad.userchanged = 1B
   ss.star[0].slope.userchanged = 1B
   ss.telescope[*].gamma.value += coeffs[0]
endif else if ss.star[0].slope.fit then begin
   coeffs = poly_fit(alltime-ss.rvepoch, allrv, 1, yfit=yfit)
   allrv -= yfit
   ss.star[0].slope.value = coeffs[1]
   ss.telescope[*].gamma.value += coeffs[0]
   ss.star[0].slope.userchanged = 1B
endif
ss.telescope[*].gamma.userchanged = 1B

if ss.ntel gt 0 then begin
   ss.planet[*].k.value = sqrt(2d0)*stddev(allrv)
   ss.planet[*].k.userchanged=1B
endif

;; read in astrometry files
if nastrom gt 0 then begin
   ss.astrom[*].astromptrs = ptrarr(nastrom,/allocate_heap)
   for i=0L, nastrom-1 do begin
      *(ss.astrom[i].astromptrs) = readastrom(astromfiles[i])
      band = (*(ss.astrom[i].astromptrs)).band
      ss.astrom[i].bandndx = where(ss.band[*].name eq band)
      ss.ndata += n_elements((*(ss.astrom[i].astromptrs)).bjd)
   endfor
endif

for i=0L, nastrom-1 do begin

   astrom = *(ss.astrom[i].astromptrs)
   if ~keyword_set(astrom.userhopa) then begin
      ss.star[0].parallax.fit = 1B
      ss.star[0].parallax.derive=1B
      ss.star[0].ra.fit = 1B
      ss.star[0].ra.derive=1B
      ss.star[0].dec.fit = 1B
      ss.star[0].dec.derive=1B
      ss.star[0].pmra.fit = 1B
      ss.star[0].pmra.derive=1B
      ss.star[0].pmdec.fit = 1B
      ss.star[0].pmdec.derive=1B
      ss.star[0].rvabs.fit = 1B
      ss.star[0].rvabs.derive=1B

      ss.astrom[i].raoffset.fit = 1B
      ss.astrom[i].raoffset.derive=1B
      ss.astrom[i].decoffset.fit = 1B
      ss.astrom[i].decoffset.derive=1B
      ;ss.astrom[i].phottobary.derive = 1B
      ss.astrom[i].astromscale.fit = 1B
      ss.astrom[i].astromscale.derive=1B
   endif
endfor

;; make the prior array
;priors = [-1,-1,-1,-1,-1]
; ************************** Make the prior array ****************************
prior0 = create_struct('name','',$  ;; the parameter name of the prior
                       'map',[-1,-1,-1,-1,-1,-1],$               ;; a map to the parameter this applies to (need not be unique)
                       'value',[-1d0,-1d0,-1d0,-1d0,-1d0,-1d0],$ ;; either a value or a map to another parameter
                       'lowerbound',-!values.d_infinity,$    ;; if value < lowerbound, loglike = -!values.d_infinity
                       'upperbound',!values.d_infinity,$     ;; if value > upperbound, loglike = -!values.d_infinity
                       'gaussian_width',!values.d_infinity,$ ;; loglike -= 0.5*((model-value)/gaussian_width)^2
                       'posterior',ptr_new())                ;; a pointer to a posterior from which to apply the prior

;; for each input prior
if ~keyword_set(silent) then printandlog, 'These are the priors that will be applied to the fit.', logname
if ~keyword_set(silent) then printandlog, 'Those with "no prior constraint" only affect the starting values of the fit:', logname
if ~keyword_set(silent) then printandlog, '', logname
if ~file_test(priorfile) then begin
   printandlog, "A priorfile must be specified. " + priorfile + " not found.", logname
   return, -1
endif

openr, lun, priorfile, /get_lun
line = ''
lineno=0

while not eof(lun) do begin

   prior = prior0
   readf, lun, line

   ;; skip commented lines
   if strpos(line,'#') eq 0 then continue
   
   ;; strip everything after comments
   entries = strsplit((strsplit(line,'#',/extract))[0],/extract)
   
   nentries = n_elements(entries)
   prior.name = entries[0]

   ;; each line must have at least a name and value
   if nentries lt 2 or nentries gt 6 then begin
      if line ne '' and ~keyword_set(silent) then $
         printandlog, 'WARNING: line ' + strtrim(lineno,2) + ' in ' + priorfile + ' is not legal syntax '+$
                      '(NAME VALUE [UNCERTAINTY] [LOWERBOUND] [UPPERBOUND] [START]) ; ignoring: ' + line, logname
      continue
   endif 

   ;; map the variable names to the structure
   prior.map = findvariable(ss,strtrim(entries[0],2),logname=logname, silent=silent, count=0)
   if prior.map[0] eq -1 then continue ;; didn't find it, skip (findvariable will print warning)

   priorvalstr = entries[1]
   ;; see if the value is a string (create a map) or a number (use the value)
   if valid_num(entries[1]) then begin
      priorval = double(entries[1])
      prior.value = [priorval,-1d0,-1d0,-1d0,-1d0,-1d0]
   endif else begin
      prior.value = findvariable(ss,strtrim(entries[1],2),logname=logname, silent=silent, count=0) 
      if prior.value[0] eq -1 then continue ;; didn't find it, skip (findvariable will print warning)
   endelse 

   ;; if the prior is a value
   if prior.map[4] ne -1 then begin
      parameter = (*ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]]).(prior.map[4])[prior.map[5]]
   endif else if prior.map[2] ne -1 then begin 
      parameter = ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]]
   endif else if prior.map[0] ne -1 then begin
      parameter = ss.(prior.map[0])[prior.map[1]]
   endif

   ;; ***** special behaviors for specific priors *****

   ;; if there's a prior on the age, fit it  
   ;; (if it's a map to another variable, we need to fit both)
   if strpos(strupcase(prior.name), 'AGE') eq 0 then begin
      parameter.fit = 1B
      if prior.value[2] ne -1 then ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].fit = 1B
   endif

   ;; if fbol prior is supplied, then we must derive it and fit the distance
   if prior.name eq 'fbol' then begin
      if ~keyword_set(silent) then printandlog, 'Fbol prior supplied, fitting distance', logname
      parameter.derive = 1B
      ss.(prior.map[0])[prior.map[1]].distance.fit = 1B
      ss.(prior.map[0])[prior.map[1]].parallax.derive = 1B
   endif
   
   ;; if the circular flag is set, don't apply priors on e/omega
   if prior.map[3] ne -1 then begin
      if circular[prior.map[3]] then begin ;; if fixed to be circular
         if prior.name eq 'e' or prior.name eq 'omega' or prior.name eq 'omegadeg' or $
            prior.name eq 'secosw' or prior.name eq 'sesinw' or prior.name eq 'esinw' or $
            prior.name eq 'ecosw' and circular[prior.map[3]] then begin
            if ~keyword_set(silent) then begin
               printandlog, "WARNING: Prior supplied on '" + $
                            prior.name + "' but planet is specified to be circular; not applying prior", logname
            endif
            continue
         endif
      endif
   endif
   ;; **************************************************
      

   ;; if the value is a map to a variable name, then get its value from that variable
   if prior.value[5] ne -1 then begin
      priorval = (*ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]]).(prior.value[4])[prior.value[5]].value
   endif else if prior.value[3] ne -1 then begin
      priorval = ss.(prior.value[0])[prior.value[1]].(prior.value[2])[prior.value[3]].value
   endif else if prior.value[1] ne -1 then begin
      priorval = ss.(prior.value[0])[prior.value[1]].value
   endif else priorval = prior.value[0]

   if n_elements(priorval) ge 2 then begin
      if prior.gaussian_width eq 0 then parameter.value = priorval
      priorval = median(priorval)
   endif

   if nentries ge 3 then begin
      if double(entries[2]) eq 0d0 then begin
         if ~parameter.fit then begin
            ;; we can't fix parameters that are not fit
            if ~keyword_set(silent) then printandlog, prior.name + ' is not fit. Only fitted parameters can be fixed. Ignoring constraint', logname
         endif else begin
            parameter.fit = 0B
            if ~keyword_set(silent) then printandlog, prior.name + ' = ' + priorvalstr + ' (fixed)', logname
            prior.gaussian_width = 0d0
         endelse
      endif else begin
         prior.gaussian_width = double(entries[2]) 
         ss.ndata++
      endelse
   endif else prior.gaussian_width = -1

   ;; if a lower bound is specified, use it
   if nentries ge 4 then begin
      if strupcase(entries[3]) eq '-INF' then prior.lowerbound = -!values.d_infinity $
      else prior.lowerbound = double(entries[3])
   endif
   
   ;; if an upper bound is specified, use it
   if nentries ge 5 then begin
      if strupcase(entries[4]) eq 'INF' then prior.upperbound = !values.d_infinity $
      else prior.upperbound = double(entries[4])
   endif

   ;; if a starting value is specified
   ;; (different from the prior value), use it
   if nentries eq 6 then begin
      startval = double(entries[5])
   endif else startval = priorval

   ;; output the priors applied to the log/screen
   if finite(prior.lowerbound) or finite(prior.upperbound) then begin
      boundtxt = '; bounded between ' + strtrim(prior.lowerbound,2) + ' and ' +  strtrim(prior.upperbound,2)
   endif else boundtxt = ''

   if ~keyword_set(silent) then begin
      if prior.gaussian_width lt 0d0 then begin
         printandlog, prior.name + ' = ' + strtrim(startval,2) + ' (no prior constraint)' + boundtxt,logname
      endif else if prior.gaussian_width gt 0 then begin
         printandlog, prior.name + ' = ' + priorvalstr + ' +/- ' + strtrim(prior.gaussian_width,2) + boundtxt, logname
      endif
      if startval ne priorval then begin
         printandlog, 'NOTE: ' + prior.name + ' starts at ' + strtrim(startval,2), logname
      endif
   endif

   if parameter.userchanged and (finite(parameter.priorwidth) or finite(parameter.lowerbound) or finite(parameter.upperbound)) then begin
      printandlog, 'WARNING: ' + prior.name + ' already set. Overwriting', logname
   end

   ;; flag that this parameter has been changed by the user
   ;; so it doesn't get overwritten later
   parameter.userchanged = 1B
   parameter.value = startval
   parameter.prior = priorval
   parameter.priorwidth = prior.gaussian_width
   parameter.lowerbound = prior.lowerbound
   parameter.upperbound = prior.upperbound

   if n_elements(*parameter.prior_new) eq 0 then begin
      parameter.prior_new = ptr_new(prior,/allocate_heap)
   endif else begin
      *parameter.prior_new = [(*(parameter.prior_new)),prior]
   endelse

   ;; warn the user that they may be specifying a useless prior, 
   ;; which is likely unintended
   if ~parameter.fit and ~parameter.derive and ~keyword_set(silent) then begin
      printandlog, "WARNING: Prior supplied on '" + $
                   prior.name + "' but it is neither fitted or derived. Not applying prior, " + $
                   'but it will be used to derive the starting parameters if possible.', logname
   endif

   ;; assign the parameter back to the structure
   if prior.map[4] ne -1 then begin
      (*ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]]).(prior.map[4])[prior.map[5]] = parameter
   endif else if prior.map[2] ne -1 then begin 
      ss.(prior.map[0])[prior.map[1]].(prior.map[2])[prior.map[3]] = parameter
   endif else if prior.map[0] ne -1 then begin
      ss.(prior.map[0])[prior.map[1]] = parameter
   endif

   ;; user only wants to change the starting value
   ;; don't create a prior entry for it
   if nentries eq 2 then continue

   ;; add it to the prior array
   if n_elements(priors) eq 0 then priors = prior $
   else priors = [priors,prior]

endwhile
if n_elements(priors) gt 0 then (*ss.priors) = priors

;; was a prior on Tc and period given? If not, use BLS (if
;; ntransits>1) or Lomb-Scargle (on RVs)
;; aspirational stub -- not implemented/tested yet
noephemeris = 0
if noephemeris then begin
   printandlog, "WARNING: no ephemeris (tc/period) given in the priors.",logname
   printandlog, "Attempting to determine ephemeris, but EXOFAST is not designed to search for planets --", logname
   printandlog, "especially for multi-planet systems.", logname
   printandlog, "Supplying an accurate ephemeris in the prior file is STRONGLY recommended.", logname

   G = constants.GMSun/constants.RSun^3*constants.day^2 ;; R_sun^3/(m_sun*day^2)
   minperiod = sqrt(4d0*!dpi^2*ss.star.rstar.value^3/(G*ss.star.mstar.value)) ;; corresponds to a/rstar=1 
   maxperiod = dataspan ;; can't see periodicities larger

   printandlog, "Searching for periods between " + strtrim(minperiod,2) + " and " + strtrim(maxperiod,2) + " days", logname
   if ntran gt 1 then begin
      printandlog, "Using BLS to determine ephemeris", logname
      bestperiod = exofast_bls(time,flux,minperiod=minperiod,maxperiod=maxperiod, tc=tt,periods=periods,ndx=ndx)
      ss.planet[i].period.value = periods[ndx]
      ss.planet[i].tt.value = tt[ndx]
   endif else if ntel gt 0 then begin
      printandlog, "Using Lomb-Scargle to determine ephemeris", logname
      bestperiod = exofast_lombscargle(time,rv,err,minperiod=minperiod, maxperiod=maxperiod, $
                                       bestpars=bestpars,slope=ss.star[0].slope.fit,quad=ss.star[0].quad.fit)
      ss.planet[i].period.value = bestpars[1]
      ss.planet[i].tc.value = bestpars[0]
   endif
endif

changed = where(ss.planet.tc.userchanged)
if changed[0] ne -1 then begin
   minbjd = min(ss.planet[changed].tc.prior,maxbjd)
   if minbjd lt minallbjd then minallbjd = minbjd
   if maxbjd gt maxallbjd then maxallbjd = maxbjd
endif
changed = where(ss.planet.tp.userchanged)
if changed[0] ne -1 then begin
   minbjd = min(ss.planet[changed].tp.prior,maxbjd)
   if minbjd lt minallbjd then minallbjd = minbjd
   if maxbjd gt maxallbjd then maxallbjd = maxbjd
endif
changed = where(ss.planet.tt.userchanged)
if changed[0] ne -1 then begin
   minbjd = min(ss.planet[changed].tt.prior,maxbjd)
   if minbjd lt minallbjd then minallbjd = minbjd
   if maxbjd gt maxallbjd then maxallbjd = maxbjd
endif

dataspan = (maxallbjd - minallbjd)
if dataspan gt 1d5 then begin
   printandlog, '', logname
   printandlog, "WARNING: data/priors span " + strtrim(dataspan/365.25d0,2) + $
                " years. Make sure all data have a consistent epoch (e.g., you're not using a mix of MJD and JD).", logname
   printandlog, "type '.con' to ignore and continue", logname
endif
free_lun, lun

if ~keyword_set(silent) then printandlog, '', logname

;; do we have enough information to derive the distance?
;if (where(priorname eq 'distance'))[0] ne -1 
;priors = priors[*,1:*]
*(ss.priors) = priors

;; determine the epoch for each observation
for i=0, ntran-1 do begin

;   if i eq 0 then begin
;      ;; variations are defined relative to the first transit
;      if keyword_set(tivs) then begin
;         ss.transit[i].tiv.fit = 0
;         ss.transit[i].tiv.derive = 0
;      endif
;      if keyword_set(tdeltavs) then begin
;         ss.transit[i].tdeltav.fit = 0
;         ss.transit[i].tdeltav.derive = 0
;      endif
;   endif
   
   for j=0, nplanets-1 do begin

      

;      ss.transit[i].epoch = min(round((mean((*(ss.transit[i].transitptrs)).bjd) - ss.planet[ss.transit[i].pndx].tc.value)/ss.planet[ss.transit[i].pndx].period.value))

      if n_elements(nvalues) ne 0 then begin
         tc = median(ss.planet[j].tc.value)
         period = median(ss.planet[j].period.value)
         if period eq 0 then period = median(10^ss.planet[j].logp.value)
      endif else begin
         tc = ss.planet[j].tc.value
         period = ss.planet[j].period.value
         if period eq 0 then period = 10^ss.planet[j].logp.value
      endelse
      
      minbjd = min((*(ss.transit[i].transitptrs)).bjd, max=maxbjd)
      

      ss.epochs[i,j] = round(((maxbjd+minbjd)/2d0 - tc)/period)

      
;      t14 = ??
;      t14s = ??

;      thistc = tc + epoch*period
;      thists = 
;      if (minbjd le thistc + t14 and maxbjd ge thistc - t14) or $
;         (minbjd le thists + t14s and maxbjd ge thists - t14s) then $
;            ;; the ith transit covers the jth planet; add it to the bitmask
;            transit[i].bitmask += 2ULL^j 
;

      ss.transit[i].epoch[j] = ss.epochs[i,j]

   endfor
endfor

;; creates an array of indicies into the stellar structure to map which parameters should be fit
;; fit[*,0] indexes the object [star=0, planet=1, band=2, telescope=3, or transit=4]
;; fit[*,1] each object can have any number of copies. 
;; This indexes which copy (e.g., planet b=0, planet c=1 or B band=0, V band=1) 
;; fit[*,2] indexes the parameter of the object (e.g., Teff=0, [Fe/H]=1)
;; this assumes a certain structure of the parameters... is that ok?
tofit = [-1,-1,-1,-1,-1]
for i=0L, n_tags(ss)-1 do begin
   for j=0L, n_elements(ss.(i))-1 do begin
      for k=0L, n_tags(ss.(i)[j])-1 do begin

         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin
            if ptr_valid(ss.(i)[j].(k)) then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'fit') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].fit then tofit = [[tofit],[i,j,k,l,m]]
                           if (*(ss.(i)[j].(k))).(l)[m].fit or (*(ss.(i)[j].(k))).(l)[m].derive then ss.npars++
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            ;; and this captures everything else
            if tag_exist(ss.(i)[j].(k),'fit') then begin
               if ss.(i)[j].(k).fit then tofit = [[tofit],[i,j,k,-1,-1]]
               if ss.(i)[j].(k).fit or ss.(i)[j].(k).derive then ss.npars++
            endif
         endif
      endfor
   endfor
endfor
tofit = tofit[*,1:*]
*(ss.tofit) = tofit
ss.nchains = n_elements((*ss.tofit)[0,*])*2L

;; populate the best values
if arg_present(best) then pars2str, best, ss, /best

;; don't do these when creating the MCMC structure
if n_elements(ss.star[0].mstar.value) eq 1 then begin

   ;; derive all step parameters
   if not pars2step(ss) then begin
;;      printandlog, 'ERROR: The isochrones are not applicable here; refine your priors.', logname
      return, 0
   endif

   ;; use this to require planets to stay in the same order as they start
   if ss.nplanets gt 0 then begin
      priorsupplied = where(ss.planet.period.value ne 0d0)
      if priorsupplied[0] ne -1 then ss.planet[priorsupplied].logp.value = alog10(ss.planet[priorsupplied].period.value)
      ss.planetorder = sort(ss.planet.logp.value)
   endif

   ;; return an error if the starting stellar system is not allowed
   if step2pars(ss,/verbose,/changedefaults) eq -1 then begin
      printandlog, 'ERROR: starting values for stellar system not allowed; refine priors', logname
      return, -1
   endif

   ;; calculate the right period stepping scale for AMOEBA based on the
   ;; input data
   perscale, ss
endif else begin
   ;; use this to require planets to stay in the same order as they start
   if ss.nplanets gt 0 then ss.planetorder = sort(ss.planet.period.value[0])
endelse

;; load stellar structure into common block
;; must do it here because threads can't pass structures
if n_elements(chi2func) eq 1 then junk = call_function(chi2func, /loadss, ss0=ss)

if (ss.mistsedfile ne '' or ss.fluxfile ne '') and $
   (nastrom eq 0 and ~finite(ss.star[0].parallax.priorwidth)) then  begin
   printandlog, 'WARNING: Fitting an SED without providing parallax information will essentially determine a photometric parallax.', logname
endif

ss.ndata += 3*total(ss.mist)
ss.ndata += 3*total(ss.yy)
ss.ndata += 3*total(ss.parsec)
ss.ndata += 2*total(ss.torres)

return, ss

end
