;; replace environment variables with their values
function expand_path, string

  directories = strsplit(string,path_sep(),/extract)
  for i=0L, n_elements(directories)-1 do begin
     environment_variable_start = strpos(string,'$')
     if environment_variable_start eq 0 then begin
        envvar = strmid(string,1) ;; trim the $
        directories[i] = getenv(envvar)
     endif
  endfor
  return, strjoin(directories,path_sep())

end
