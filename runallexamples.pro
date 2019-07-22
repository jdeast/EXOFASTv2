;; runs short versions of all the examples (mostly as a very rough
;; unit test to make sure I didn't break anything).

pro runallexamples, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps = 100
if n_elements(nthin) eq 0 then nthin = 1

fithat3, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fithat3notran, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt1, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt4, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt6, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt12, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt16, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitwasp76, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fit8361, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitgj9827, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbos
fithd106315, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fithat3_comparestar, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose

;; test license free uses
spawn, 'gdl -e "fithat3, maxsteps=' + strtrim(maxsteps,2) + '"'

argfile = filepath('hat3.args',$
                   root_dir=getenv('EXOFAST_PATH'),$
                   subdir=['examples','hat3_nolicense'])
vmfile = filepath('exofastv2.sav',$
                   root_dir=getenv('EXOFAST_PATH'))

spawn, 'idl -arg ' + argfile + ' -vm=' + vmfile

stop
end
