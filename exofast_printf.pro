;; this is a replacement for printf that uses a (slow)
;; workaround if a demo license is being used
;; does not properly escape double quotes inside str

pro exofast_printf, lun, str, filename=filename, new=new

  if lmgr(/demo) then begin
     if n_elements(filename) eq 0 then $
        message, 'filename must be specified in DEMO mode'
     if keyword_set(new) then pipetxt = ' > ' $
     else pipetxt = ' >> '
     spawn, 'echo "' + str + '"' + pipetxt + filename
  endif else begin
     if keyword_set(new) then begin 
        if n_elements(filename) eq 0 then $
           message, 'filename must be specified when /NEW is set'     
        openw, lun, filename, /get_lun
     endif
     printf, lun, str
     if keyword_set(close) then free_lun, lun
  endelse

end
