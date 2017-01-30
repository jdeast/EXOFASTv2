function readtran, filename

if not file_test(filename) then message, 'Transit file (' + filename + ') does not exist'

if n_elements(strsplit(filename,'.',/extract)) lt 3 then message, 'filename (' + filename + ') must have format nYYYYMMDD.FILTER.TELESCOPE.whateveryouwant'

;; Read the transit data file into a structure
;; (with an arbitary number of detrending variables)
band = (strsplit(filename,'.',/extract))(1)

if band eq 'u' then begin
   band = 'Sloanu'
   bandname = "u"
endif else if band eq 'g' then begin
   band = 'Sloang'
   bandname = "g"
endif else if band eq 'r' then begin
   band = 'Sloanr'
   bandname = "r"
endif else if band eq 'i' then begin
   band = 'Sloani'
   bandname = "i"
endif else if band eq 'z' then begin
   band = 'Sloanz'
   bandname = "z"
endif else bandname = band

allowedbands = ['U','B','V','R','I','J','H','K',$
                'Sloanu','Sloang','Sloanr','Sloani','Sloanz',$
                'Kepler','CoRoT','Spit36','Spit45','Spit58','Spit80',$
                'u','b','v','y']

if (where(allowedbands eq band))[0] eq -1 then message, 'Filter (' + band + ') not allowed'

line = ""
openr, lun, filename, /get_lun
readf, lun, line
entries = double(strsplit(line,/extract))
ncol = n_elements(entries)
nrow = file_lines(filename)

;; rewind file to beginning
point_lun, lun, 0
array = dblarr(ncol,nrow)
readf, lun, array
free_lun, lun

bjd = transpose(array[0,*])
flux = transpose(array[1,*])
err = transpose(array[2,*])
if ncol gt 3 then begin
   da = array[3:ncol-1,*]
   ndetrend = ncol-3
endif else begin
   da = 0d0
   ndetrend=0
endelse

night = strmid(filename,1,4)+'-'+strmid(filename,5,2)+'-'+strmid(filename,7,2)
label = (strsplit(filename,'.',/extract))(2) + ' UT ' + night + ' ('+ bandname + ')'
transit=create_struct('bjd',bjd,'flux',flux,'err',err,'band',band,'ndx',0,$
                      'epoch',0.0,'detrendadd',da,'detrendmult',da,'label',$
                      label,'ndetrend',ndetrend)
return, transit

end
