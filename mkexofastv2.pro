;; compiles exofastv2 as an IDL save file for a virtual machine
;; (license free use)
pro mkexofastv2

;; create auxialliary functions
resolve_all, resolve_procedure='mksed',/cont,/quiet
save, /routines, filename='mksed.sav'

resolve_all, resolve_procedure='mkprior',/cont,/quiet
save, 'getpriorline','mkprior','strsplit','tag_exist',/routines, filename='mkprior.sav'

;; create the main EXOFASTv2 executeable
resolve_all, resolve_procedure='exofastv2',$ ;; the main routine
             resolve_function=['exofast_chi2v2','exofast_random'],$ ;; these are run with call_function and do not get resolved automatically
             /cont,/quiet
save, /routines, filename='exofastv2.sav'

;; create the HTML help pages
mk_html_help, './', 'exofastv2.html'

end
