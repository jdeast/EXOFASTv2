pro exofast_fixps, psname

openr, lun, psname, /get_lun
line = ''
while not eof(lun) do begin
   readf, lun, line
   if strpos(line,'%%BoundingBox:') eq 0  then begin
      entries = strsplit(line,/extract)
      docxy = long(entries[1:4])
      newdocxy = docxy
   endif else if strpos(line,'%%PageBoundingBox:') eq 0  then begin
      entries = strsplit(line,/extract)
      pagexy = long(entries[1:4])
      if pagexy[3] gt docxy[3] then newdocxy[3] = pagexy[3]
   endif
endwhile
free_lun, lun

;; no change necessary
if max(abs(docxy - newdocxy)) eq 0d0 then return

newpsname = psname + '.bak.ps'
file_copy, psname, newpsname,/overwrite

;; rewrite a new file replacing the bounding box
openw, newlun, newpsname, /get_lun
openr, oldlun, psname, /get_lun
while not eof(oldlun) do begin
   readf, oldlun, line
   if strpos(line,'%%BoundingBox:') eq 0  then begin
      printf, newlun, '%%BoundingBox: ' + strjoin(strtrim(newdocxy,2),' ')
   endif else printf, newlun, line
endwhile
free_lun, newlun
free_lun, oldlun

file_move, newpsname, psname,/overwrite

end
      
