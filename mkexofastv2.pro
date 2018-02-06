;; compiles exofastv2 as an IDL save file for a virtual machine
;; (license free use)
pro mkexofastv2

resolve_all, resolve_procedure='exofastv2',$ ;; the main routine
             resolve_function=['exofast_chi2v2','exofast_random'],$ ;; these are run with call_function and do not get resolved automatically
             /cont,/quiet

save, /routines, filename='exofastv2.sav'

mk_html_help, './', 'exofastv2.html'

end
