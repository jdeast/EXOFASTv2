pro getfilter, filterid, redo=redo

constants = mkconstants()
c = constants.c*1d4 ;; microns/s

idlname = str_replace(filterid,'/','_') + '.idl'
idlfilepath = filepath(idlname,root_dir=getenv('EXOFAST_PATH'),subdir=['sed','filtercurves'])
if file_test(idlfilepath) and (~keyword_set(redo)) then return

xmlname = str_replace(filterid,'/','_') + '.xml'
xmlfilepath = filepath(xmlname,root_dir=getenv('EXOFAST_PATH'),subdir=['sed','filtercurves'])

if not file_test(xmlfilepath) then begin
   print, "no filter with ID " + filterid + ', Attempting to retrieve from svo2'
   url = 'http://svo2.cab.inta-csic.es/theory/fps3/fps.php?ID=' + filterid
   junk = wget(url, filename=xmlfilepath)

   if not file_test(xmlfilepath) then begin
      print, 'retrieval failed, skipping ' + filterid
      return
   endif
endif          

wavelength=[]
transmission=[]
line = ''
openr, lun, xmlfilepath, /get_lun
while not eof(lun) do begin
   readf, lun, line
   if strpos(line,'WavelengthEff') ne -1 then begin
      weff = double((strsplit((strsplit(line,'value="',/extract,/regex))[1],'"',/extract))[0])/1d4 ;; microns
   endif else if strpos(line,'WidthEff') ne -1 then begin
      widtheff = double((strsplit((strsplit(line,'value="',/extract,/regex))[1],'"',/extract))[0])/1d4 ;; microns
   endif else if strpos(line,'name="ZeroPoint"') ne -1 then begin
      zero_point = double((strsplit((strsplit(line,'value="',/extract,/regex))[1],'"',/extract))[0]) ;; Jy
   endif else if strpos(line,'<TR>') ne -1 then begin
      readf, lun, line
      wavelength = [wavelength,double((strsplit(line,'><',/extract))[2])/1d4]
      readf, lun, line
      transmission = [transmission,double((strsplit(line,'><',/extract))[2])]
   endif
endwhile
free_lun, lun

;; wavelength scale on which to interpolate
intwave = dindgen(24000L)/1000d0 + 0.1d0

;; pad ends with zeros to help interpolation
lowwave = where(intwave lt wavelength[0])
if lowwave[0] ne -1 then begin
   transmission = [lowwave*0d0,transmission]
   wavelength = [intwave[lowwave],wavelength]
endif
highwave = where(intwave gt wavelength[n_elements(wavelength)-1])
if highwave[0] ne -1 then begin
   transmission = [transmission,highwave*0]
   wavelength = [wavelength,intwave[highwave]]
endif

inttran = interpol(transmission,wavelength,intwave)

filter = create_struct('id',filterid,$
                       'weff',weff,$ ;; effective wavelength in microns
                       'widtheff',widtheff,$              ;; effective width in microns
                       'zero_point',1d-23*zero_point*c/weff,$ ;; zero point (Jy -> erg/s/cm^2)
                       'ab_offset',0d0,$
                       'wavelength',intwave,$ ;; wavelength in microns
                       'transmission',inttran) ;; transmission/wavelength in fraction
;                       'transmission',inttran/(1000d0)) ;; transmission/wavelength in fraction

;; flux from magnitude
;flux = filter.zero_point*10^(-0.4d0*(mag + filter.ab_offset))

;; flux from sed
;nwaves = 24000
;wavelength = 
;transmission = interpol([0d0,filter.transmission,0d0],[0.1d0,filter.wavelength/1d4,24.099d0],wavelength)
;flux = total(sed*transmission)/filter.widtheff*(wavelength[1]-wavelength[0])

save, filter, filename=idlfilepath
print, 'successfully retrieved ' + filterid

end

pro getallfilters

filterids = ['GAIA/GAIA2r.Gbp',$
             'GAIA/GAIA2r.Grp',$
             'GAIA/GAIA2r.G',$
             'GAIA/GAIA3.Gbp',$
             'GAIA/GAIA3.Grp',$
             'GAIA/GAIA3.G',$
             'Keck/NIRC2.Ks',$
             'Keck/NIRC2.Brgamma',$
             '2MASS/2MASS.J',$
             '2MASS/2MASS.H',$
             '2MASS/2MASS.Ks',$
             'WISE/WISE.W1',$
             'WISE/WISE.W2',$
             'WISE/WISE.W3',$
             'WISE/WISE.W4',$
             'TESS/TESS.Red',$
             'Kepler/Kepler.K',$
             'WIYN/NESSI.NB467',$ ;; Commonly used AO bands to split close binaries 
             'WIYN/NESSI.NB562',$ ;; Commonly used AO bands to split close binaries 
             'WIYN/NESSI.NB716',$ ;; Commonly used AO bands to split close binaries 
             'WIYN/NESSI.NB832',$ ;; Commonly used AO bands to split close binaries 
             'Gemini/Zorro.EO_466',$ ;; Commonly used AO bands to split close binaries 
             'Gemini/Zorro.EO_562',$ ;; Commonly used AO bands to split close binaries 
             'Gemini/Zorro.EO_716',$ ;; Commonly used AO bands to split close binaries 
             'Gemini/Zorro.EO_832',$ ;; Commonly used AO bands to split close binaries 
             'WIYN/NESSI.u',$
             'WIYN/NESSI.g',$
             'WIYN/NESSI.r',$
             'WIYN/NESSI.i',$
             'WIYN/NESSI.z',$
             'Misc/APASS.B',$
             'Misc/APASS.V',$
             'Generic/Cousins.I',$
             'Misc/APASS.sdss_g',$
             'Misc/APASS.sdss_r',$
             'Misc/APASS.sdss_i',$
             'Misc/APASS.sdss_z']
             
for i=0L, n_elements(filterids)-1 do begin
   getfilter, filterids[i],/redo
endfor

end
