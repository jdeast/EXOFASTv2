function readrv_detrend, detrendpar=detrend, filename

label = (strsplit(filename,'.',/extract))(1)

line = ''
openr, lun, filename,  /get_lun
readf, lun, line
if strpos(line,'SB2_') ne -1 then begin
   planet = long((strsplit(line,'SB2_',/extract,/regex))[1])
   sb2line = 1B
   readf, lun, line
endif else begin
   planet = -1L
   sb2line = 0B
endelse

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
   ;; now read the first data line
   readf, lun, line 
endif else begin
   header = 0
   mult = [-1]
   nmult=0
endelse

entries = double(strsplit(line,/extract))
ncol = n_elements(entries)
;; if no header, assume all columns are additive
if not header then begin
   if ncol le 3 then begin
      add = [-1]
      nadd = 0L
   endif else begin
      nadd = ncol-3
      add = [-1,3L + lindgen(nadd)]
   endelse
endif

if ncol lt 3 then begin
   message, 'RV file (' + filename + ') must contain at least 3 white-space delimited columns (BJD_TDB RV RVerr).',/continue
   message, 'If fitting an SB2, the first line must be "SB2_#, where "#" is the index of the corresponding planet',/continue
   message, 'The first line (or second line if fitting an SB2) can be a header, if it starts with "#"'
   message, 'The first line is ' + line
endif

if nadd+nmult+3 ne ncol then begin
   message, 'Mismatch between header line and data lines for ' + filename + '. The header MUST have one white-space delimited entry per column.'
endif

nrow = file_lines(filename) - header - sb2line

;; rewind file to beginning
point_lun, lun, 0
if sb2line then readf, lun, line
if header then readf, lun, line
array = dblarr(ncol,nrow)
readf, lun, array
free_lun, lun

bjd = transpose(array[0,*])
rv = transpose(array[1,*])
err = transpose(array[2,*])

bi = err + !values.d_nan
bierr = err  + !values.d_nan


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
detrendaddpars.label = 'RVC' + strtrim(indgen(nadd > 1),2)
detrendaddpars.latex = 'RVC_{' + strtrim(indgen(nadd > 1),2) + '}'
if nadd eq 0 then begin
   detrendaddpars.fit = 0
   detrendaddpars.derive = 0
endif

detrendmultpar = detrend
detrendmultpar.description = 'Multiplicative detrending coeff'
detrendmultpars = replicate(detrendmultpar,nmult > 1)
detrendmultpars.label = 'RVM' + strtrim(indgen(nmult > 1),2)
detrendmultpars.latex = 'RVM_{' + strtrim(indgen(nmult > 1),2) + '}'
if nmult eq 0 then begin
   detrendmultpars.fit = 0
   detrendmultpars.derive = 0
endif

return, create_struct('bjd',bjd,'rv',rv,'err',err, 'bi',bi,'bierr',bierr,$
                      'label',label,'residuals',dblarr(n_elements(bjd)),$
                      'rm',dblarr(n_elements(bjd)), 'planet',planet,$
                      'nadd',nadd,'nmult',nmult,$
                      'detrendadd',da, 'detrendmult',dm,$
                      'detrendaddpars',detrendaddpars, 'detrendmultpars',detrendmultpars)

end
