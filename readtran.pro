function readtran, filename, detrendpar=detrend, nplanets=nplanets

if n_elements(nplanets) eq 0 then nplanets = 1d0

if not file_test(filename) then message, 'Transit file (' + filename + ') does not exist'

basename = file_basename(filename)
if n_elements(strsplit(basename,'.',/extract)) lt 3 then message, 'filename (' + basename + ') must have format nYYYYMMDD.FILTER.TELESCOPE.whateveryouwant (see readtran.pro for details)'

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

line = ""
openr, lun, filename, /get_lun
readf, lun, line

breakptdates = [-1]
breakptline = 0

if strpos(line,'#') eq 0 then begin
   mult = [-1]
   add = [-1]
   header = 1
   nadd = 0L
   nmult=0L
   entries = strsplit(line,'# '+string(9B),/extract)
   ncol = n_elements(entries)
   for i=3, ncol-1 do begin
      if strpos(entries[i],'M') eq 0 then begin
         ;; if it's preceeded by an M, make it a multiplicative detrending
         mult = [mult,i]
         nmult+=1L
      endif else begin
         ;; otherwise, make it additive
         add = [add,i]
         nadd+=1L
      endelse
   endfor

   ;; now read the next line
   readf, lun, line 
   if strpos(line,'#') eq 0 then begin
      breakptline = 1
      ;; then this is a special line denoting breakpoints for the spline fit of the OOT lightcurve
      breakptdates = double(strsplit(line,'# ',/extract))
      readf, lun, line ;; now read the first data line
   endif 
endif else begin
   header = 0
   mult = [-1]
   nmult=0
endelse

;; for TESS lightcurves, use the momentum dumps for break points by default
;if band eq 'TESS' and breakptdates[0] eq -1 then begin
;   readcol, getenv('EXOFAST_PATH') + path_sep() + 'tess_momentumdumps.txt', breakptdates, format='d'
;endif

entries = double(strsplit(line,/extract))
ncol = n_elements(entries)
;; if no header, assume all additive
if not header then begin
   if ncol le 3 then begin
      add = [-1]
      nadd = 0L
   endif else begin
      nadd = ncol-3
      add = [-1,3L + lindgen(nadd)]
   endelse
endif

if ncol lt 3 then message, 'Transit file (' + filename + ') must contain at least 3 white-space delimited columns (BJD_TDB flux err). Comments are not allowed. The first line is ' + line

nrow = file_lines(filename) - header - breakptline

;; rewind file to beginning
point_lun, lun, 0
if header then readf, lun, line
if breakptline then readf, lun, line
array = dblarr(ncol,nrow)
readf, lun, array
free_lun, lun

bjd = transpose(array[0,*])
flux = transpose(array[1,*])
err = transpose(array[2,*])

;; translate break point times to break point indices
breakpts = [-1]
if breakptdates[0] ne -1 then begin
   for i=0L, n_elements(breakptdates)-1 do begin
      junk = min(abs(bjd-breakptdates[i]),ndx)
      if ndx ne 0 and ndx ne n_elements(bjd)-1 then breakpts = [breakpts,ndx]
   endfor
endif

;; trim out the initial -1
if n_elements(breakpts) gt 1 then begin
   breakpts = breakpts[1:n_elements(breakpts)-1]

   ;; remove any duplicates, arrange indices in order
   sorted = sort(breakpts)
   breakpts = breakpts[sorted]
   breakpts = breakpts[uniq(breakpts)]

   ;; remove adjacent break points (which breaks the spline)
   good = where(breakpts-shift(breakpts,1) ne 1)
   breakpts = breakpts[good]
endif

;; zero average the detrending parameters
array -= transpose(total(array,2)/n_elements(array[0,*])##replicate(1d0,n_elements(array[0,*])))

;; rescale the detrending parameters to the range -1 to 1
array /= transpose(max(abs(array),dim=2)##replicate(1d0,n_elements(array[0,*])))

if nadd gt 0 then da = array[add[1:nadd],*] $
else da = 0d0

if nmult gt 0 then dm = array[mult[1:nmult],*] $
else dm = 0d0

detrendaddpar = detrend
detrendaddpar.description = 'Additive detrending coeff'
detrendaddpars = replicate(detrendaddpar,nadd > 1)
detrendaddpars.label = 'C' + strtrim(indgen(nadd > 1),2)
detrendaddpars.latex = 'C_{' + strtrim(indgen(nadd > 1),2) + '}'
if nadd eq 0 then begin
   detrendaddpars.fit = 0
   detrendaddpars.derive = 0
endif

detrendmultpar = detrend
detrendmultpar.description = 'Multiplicative detrending coeff'
detrendmultpars = replicate(detrendmultpar,nmult > 1)
detrendmultpars.label = 'M' + strtrim(indgen(nmult > 1),2)
detrendmultpars.latex = 'M_{' + strtrim(indgen(nmult > 1),2) + '}'
if nmult eq 0 then begin
   detrendmultpars.fit = 0
   detrendmultpars.derive = 0
endif

residuals = flux*0d0
model = flux*0d0
;model = (flux*0d0)#replicate(1d0,nplanets+1)

span = (max(bjd)+1d0 - (min(bjd)-1d0))
npretty = span*1440d0/5d0
prettytime = min(bjd) - 1d0 + dindgen(npretty)/(npretty-1)*span
;prettytime = replicate(min(bjd) - 1d0 + dindgen(npretty)/(npretty-1)*span,nplanets+1)
prettymodel = prettytime*0d0

night = strmid(basename,1,4)+'-'+strmid(basename,5,2)+'-'+strmid(basename,7,2)
label = (strsplit(basename,'.',/extract))(2) + ' UT ' + night + ' ('+ bandname + ')'

transit=create_struct('bjd',bjd,'flux',flux,'err',err,'band',band,'ndx',0,$
                      'epoch',0.0,'detrendadd',da,'detrendmult',dm,'label',$
                      label,'nadd',nadd,'nmult',nmult,$
                      'residuals',residuals, 'model',model, $
                      'prettytime',prettytime, 'prettymodel',prettymodel,$
                      'detrendaddpars',detrendaddpars, 'detrendmultpars',detrendmultpars, $
                      'breakpts',breakpts)

return, transit

end
