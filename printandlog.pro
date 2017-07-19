pro printandlog, string, logname, format=format

if n_elements(logname) ne 0 then begin
   if logname ne '' then begin
      openw, lun, logname,/append, /get_lun
      printf, lun, string, format=format
      free_lun, lun
   endif
endif

print, string, format=format

end
