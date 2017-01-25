pro fixldfiles

filters = ['u','b','v','y']
type = ['.linear.','.','.nonlin.']

for i=0, n_elements(filters)-1 do begin
   for k=0, n_elements(type)-1 do begin
      files = file_search('*.' + filters[i] + type + 'sav',count=nfiles)
      for j=0, nfiles-1 do begin
         split = strsplit(files[j],'.',/extract)
         split[3] = 'Strom' + filters[i]
         newfilename = strjoin(split,'.')
         print, 'Moving ' + files[j] + ' to ' +  newfilename
         file_move, files[j], newfilename
      endfor
   endfor
endfor   

end
