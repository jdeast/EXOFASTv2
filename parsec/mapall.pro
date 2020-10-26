pro mapall,overwrite=overwrite

files = file_search('Z*/*.DAT', count=nfiles)

;; just do the ones required for the Sun
;files = ['Z0.014Y0.273/Z0.014Y0.273OUTA1.74_F7_M001.000.DAT',$
;         'Z0.017Y0.279/Z0.017Y0.279OUTA1.74_F7_M001.000.DAT',$
;         'Z0.014Y0.273/Z0.014Y0.273OUTA1.74_F7_M001.050.DAT',$
;         'Z0.017Y0.279/Z0.017Y0.279OUTA1.74_F7_M001.050.DAT',$
;         'Z0.02Y0.284/Z0.02Y0.284OUTA1.74_F7_M001.000.DAT',$
;         'Z0.02Y0.284/Z0.02Y0.284OUTA1.74_F7_M001.050.DAT']
nfiles = n_elements(files)
shuffle = sort(randomu(seed,nfiles))
files = files[shuffle]

for i=0L, nfiles-1 do begin
   if strpos(files[i],".HB.") ne -1 then continue
   if strpos(files[i],".ADD.") ne -1 then continue
   print, 'starting ' + files[i] + ' (' + strtrim(i+1,2) + ' of ' + strtrim(nfiles,2) + ')'
   map2eep, files[i],overwrite=overwrite
endfor

end
