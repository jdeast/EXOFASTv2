;; runs short versions of all the examples (mostly as a very rough
;; unit test to make sure I didn't break anything).

pro runallexamples, debug=debug, verbose=verbose

fithat3, maxsteps=100, nthin=1, debug=debug, verbose=verbose
fitkelt4, maxsteps=100, nthin=1, debug=debug, verbose=verbose
fitgj9827, maxsteps=100, nthin=1, debug=debug, verbose=verbos
fithd106315, maxsteps=100, nthin=1, debug=debug, verbose=verbose
fitkelt6, maxsteps=100, nthin=1, debug=debug, verbose=verbose

stop
end
