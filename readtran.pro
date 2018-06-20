function readtran, filename, detrendpar=detrend

if not file_test(filename) then message, 'Transit file (' + filename + ') does not exist'

basename = file_basename(filename)
if n_elements(strsplit(basename,'.',/extract)) lt 3 then message, 'filename (' + basename + ') must have format nYYYYMMDD.FILTER.TELESCOPE.whateveryouwant (see readtran.pro for details)'

;; Read the transit data file into a structure
;; (with an arbitary number of detrending variables)
band = (strsplit(basename,'.',/extract))(1)

;; this means stromgren u can't be used -- remove
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

if strpos(line,'#') eq 0 then begin
   mult = []
   add = []
   header = 1
   entries = strsplit(line,'# '+string(9B),/extract)
   ncol = n_elements(entries)
   for i=3, ncol-1 do begin
      if strpos(entries[i],'M') eq 0 then begin
         ;; if it's preceeded by an M, make it a multiplicative detrending
         mult = [mult,i]
      endif else begin
         ;; otherwise, make it additive
         add = [add,i]
      endelse
   endfor
   readf, lun, line ;; now read the first data line
endif else begin
   header = 0
   mult = []
endelse

entries = double(strsplit(line,/extract))
ncol = n_elements(entries)
;; if no header, assume all additive
if not header then begin
   if ncol le 3 then add = [] $
   else add = 3L + lindgen(ncol-3)
endif

if ncol lt 3 then message, 'Transit file (' + filename + ') must contain at least 3 white-space delimited columns (BJD_TDB flux err). Comments are not allowed. The first line is ' + line

nadd = n_elements(add)
nmult = n_elements(mult)

nrow = file_lines(filename) - header

;; rewind file to beginning
point_lun, lun, 0
if header then readf, lun, line
array = dblarr(ncol,nrow)
readf, lun, array
free_lun, lun

bjd = transpose(array[0,*])
flux = transpose(array[1,*])
err = transpose(array[2,*])

;; zero average the detrending parameters
array -= transpose(total(array,2)/n_elements(array[0,*])##replicate(1d0,n_elements(array[0,*])))

;; rescale the detrending parameters to the range -1 to 1
array /= transpose(max(abs(array),dim=2)##replicate(1d0,n_elements(array[0,*])))

if nadd gt 0 then da = array[add,*] $
else da = 0d0

if nmult gt 0 then dm = array[mult,*] $
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

night = strmid(basename,1,4)+'-'+strmid(basename,5,2)+'-'+strmid(basename,7,2)
label = (strsplit(basename,'.',/extract))(2) + ' UT ' + night + ' ('+ bandname + ')'

transit=create_struct('bjd',bjd,'flux',flux,'err',err,'band',band,'ndx',0,$
                      'epoch',0.0,'detrendadd',da,'detrendmult',dm,'label',$
                      label,'nadd',nadd,'nmult',nmult,$
                      'residuals',residuals,'model',model, $
                      'detrendaddpars',detrendaddpars, 'detrendmultpars',detrendmultpars)


return, transit

end
