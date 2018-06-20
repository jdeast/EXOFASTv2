;; runs short versions of all the examples (mostly as a very rough
;; unit test to make sure I didn't break anything).

pro runallexamples, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(maxsteps) eq 0 then maxsteps = 100
if n_elements(nthin) eq 0 then nthin = 1

fithat3, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt4, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitgj9827, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbos
fithd106315, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fitkelt6, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose
fithat3_comparestar, maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose

stop
end
