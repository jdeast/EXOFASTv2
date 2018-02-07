;; reads an EEP file into a structure with labels from the header
function readmisteep, filename

line = ''
openr, lun, filename, /get_lun
entries = []
while not eof(lun) do begin
   readf, lun, line
   if strpos(line,'#') ne 0 then begin
      if n_elements(entries) eq 0 then entries = double(strsplit(line,/extract))$
      else entries = [[entries], [double(strsplit(line,/extract))]]
   endif else begin
      hdr = strtrim(strsplit(line,' #',/extract),2)
   endelse
endwhile
free_lun, lun

strcmd = 'eep = create_struct('
for i=0, n_elements(hdr)-1 do begin
   strcmd += '"' + hdr[i]+'",transpose(entries[' + strtrim(i,2) + ',*])'
   if i ne n_elements(hdr)-1 then strcmd += ','
endfor
strcmd += ')'
void = execute(strcmd)

return, eep

end



