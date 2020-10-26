pro fithat3_threads, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])
nthin=1

maxsteps=100
nthin=1

delay = [0,2.5d6,2.5d7]
      

;; Fit using the Torres relation
for j=0L, 2 do begin
   for i=1L, 8 do begin

      exofastv2, nplanets=1, tranpath=path+'n20070428.Sloani.KepCam.dat', $
                 rvpath=path+'HAT-3b.*.rv',priorfile=path+'hat3.torres.priors', $
                 prefix=path+'fitresults' + path_sep()+'HAT-3b.Threaded.MIST.' + string(j,i,format='(i03,".",i03)') + '.', mistsedfile='hat3.sed',$
                 maxsteps=maxsteps, $
                 nthin=nthin, circular=[1], fitrv=[1],fittran=[1], $
                 debug=debug, verbose=verbose,  nthreads=i, delay=delay[j]
   endfor
endfor

end
