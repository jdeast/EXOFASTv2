;; intended as a drop in replacement for forprint, but compatible with
;; virtual machines (no EXECUTE statement). Not all features
;; are implemented (but added a length check).

pro exofast_forprint, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, $
      v15,v16,v17,v18,TEXTOUT = textout, FORMAT = format, SILENT = SILENT, $ 
      STARTLINE = startline, NUMLINE = numline, COMMENT = comment, $
      SUBSET = subset, NoCOMMENT=Nocomment,STDOUT=stdout, WIDTH=width

if n_elements(textout) eq 0 then lun = -1 $
else if ~lmgr(/demo) then openw, lun, textout, /get_lun

nv1 = n_elements(v1)
if n_elements(v2) ne 0 and n_elements(v2) ne nv1 then message, 'dimensions of v2 ' + strtrim(n_elements(v2),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v3) ne 0 and n_elements(v3) ne nv1 then message, 'dimensions of v3 ' + strtrim(n_elements(v3),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v4) ne 0 and n_elements(v4) ne nv1 then message, 'dimensions of v4 ' + strtrim(n_elements(v4),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v5) ne 0 and n_elements(v5) ne nv1 then message, 'dimensions of v5 ' + strtrim(n_elements(v5),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v6) ne 0 and n_elements(v6) ne nv1 then message, 'dimensions of v6 ' + strtrim(n_elements(v6),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v7) ne 0 and n_elements(v7) ne nv1 then message, 'dimensions of v7 ' + strtrim(n_elements(v7),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v8) ne 0 and n_elements(v8) ne nv1 then message, 'dimensions of v8 ' + strtrim(n_elements(v8),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v9) ne 0 and n_elements(v9) ne nv1 then message, 'dimensions of v9 ' + strtrim(n_elements(v9),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v10) ne 0 and n_elements(v10) ne nv1 then message, 'dimensions of v10 ' + strtrim(n_elements(v10),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v11) ne 0 and n_elements(v11) ne nv1 then message, 'dimensions of v11 ' + strtrim(n_elements(v11),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v12) ne 0 and n_elements(v12) ne nv1 then message, 'dimensions of v12 ' + strtrim(n_elements(v12),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v13) ne 0 and n_elements(v13) ne nv1 then message, 'dimensions of v13 ' + strtrim(n_elements(v13),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v14) ne 0 and n_elements(v14) ne nv1 then message, 'dimensions of v14 ' + strtrim(n_elements(v14),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v15) ne 0 and n_elements(v15) ne nv1 then message, 'dimensions of v15 ' + strtrim(n_elements(v15),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v16) ne 0 and n_elements(v16) ne nv1 then message, 'dimensions of v16 ' + strtrim(n_elements(v16),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v17) ne 0 and n_elements(v17) ne nv1 then message, 'dimensions of v17 ' + strtrim(n_elements(v17),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'
if n_elements(v18) ne 0 and n_elements(v18) ne nv1 then message, 'dimensions of v18 ' + strtrim(n_elements(v18),2) + ' does not match v1 (' + strtrim(nv1,2) + ')'

if n_elements(comment) ne 0 then begin
   if lmgr(/demo) and n_elements(textout) ne 0 then spawn, 'echo "' +  '" >> ' + textout $
   else printf, lun, comment
endif

for i=0, n_elements(v1)-1 do begin

   if n_elements(v18) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i],v18[i], format=format)
   endif else if n_elements(v17) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],v17[i], format=format)
   endif else if n_elements(v16) ne 0  then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i], format=format)
   endif else if n_elements(v15) ne 0  then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i],v14[i],v15[i], format=format)
   endif else if n_elements(v14) ne 0  then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i],v14[i], format=format)
   endif else if n_elements(v13) ne 0  then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i],v13[i], format=format)
   endif else if n_elements(v12) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i],v12[i], format=format)
   endif else if n_elements(v11) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i],v11[i], format=format)
   endif else if n_elements(v10) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i],$
              v10[i], format=format)
   endif else if n_elements(v9) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],v9[i], format=format)
   endif else if n_elements(v8) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i], format=format)
   endif else if n_elements(v7) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i], format=format)
   endif else if n_elements(v6) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i],v6[i], format=format)
   endif else if n_elements(v5) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i],v5[i], format=format)
   endif else if n_elements(v4) ne 0 then begin
      str = string(v1[i],v2[i],v3[i],v4[i], format=format)
   endif else if n_elements(v3) ne 0 then begin
      str = string(v1[i],v2[i],v3[i], format=format)
   endif else if n_elements(v2) ne 0 then begin
      str = string(v1[i],v2[i], format=format)
   endif else if n_elements(v1) ne 0 then begin
      str = string(v1[i], format=format)
   endif

   if lmgr(/demo) and n_elements(textout) ne 0 then spawn, 'echo "' +  '" >> ' + textout $
   else printf, lun, str

endfor

if n_elements(textout) ne 0 and ~lmgr(/demo) then free_lun, lun

end
