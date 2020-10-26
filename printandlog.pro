pro printandlog, str, logname, format=format

if n_elements(logname) ne 0 then begin
   if logname ne '' then begin
      if file_test(file_dirname(logname)) then begin
         if ~lmgr(/demo) then begin
            openw, lun, logname, /append, /get_lun
            printf, lun, str, format=format
            free_lun, lun
         endif else begin
            spawn, 'echo "' + string(str,format=format) + '" >> ' + logname
         endelse
      endif else print, 'Warning: specified log path does not exist; not creating log'
   endif
endif

print, str, format=format

end
