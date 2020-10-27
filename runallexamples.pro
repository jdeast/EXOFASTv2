;; runs short versions of all the examples (mostly as a very rough
;; unit test to make sure I didn't break anything).

pro runallexamples, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin, nthread=nthread, compareonly=compareonly, threshhold=threshhold

if n_elements(threshhold) eq 0 then threshhold=1d0

example_dir = getenv('EXOFAST_PATH') + path_sep() + 'examples' + path_sep()

;           planet name   |  procedures |             example directory          |         ref file        |         csv file
fits = [[         'HAT-3b',    'fithat3',example_dir +        'hat3' + path_sep(),         'HAT-3b.ref.csv',        'HAT-3b.median.csv'],$ ;; cannonical example (simplified)
        [        'HATS-7b',   'fithats7',example_dir +       'hats7' + path_sep(),        'HATS-7b.ref.csv',       'HATS-7b.median.csv'],$
        [        'KELT-1b',   'fitkelt1',example_dir +       'kelt1' + path_sep(),        'KELT-1b.ref.csv',       'KELT-1b.median.csv'],$ 
        [       'KELT-4Ab',   'fitkelt4',example_dir + 'kelt4rvonly' + path_sep(),       'KELT-4Ab.ref.csv',      'KELT-4Ab.median.csv'],$ ;; RV-only example
        [        'KELT-6b',   'fitkelt6',example_dir +       'kelt6' + path_sep(),        'KELT-6b.ref.csv',       'KELT-6b.median.csv'],$ ;; LC fitting with detrending
        [       'KELT-12b',  'fitkelt12',example_dir +      'kelt12' + path_sep(),       'KELT-12b.ref.csv',      'KELT-12b.median.csv'],$ ;; TTVs, RV slope
        [       'KELT-14b',  'fitkelt14',example_dir +      'kelt14' + path_sep(),       'KELT-14b.ref.csv',      'KELT-14b.median.csv'],$ ;; TTVs
        [       'KELT-16b',  'fitkelt16',example_dir +      'kelt16' + path_sep(),       'KELT-16b.ref.csv',      'KELT-16b.median.csv'],$
        [       'KELT-17b',  'fitkelt17',example_dir +      'kelt17' + path_sep(),       'KELT-17b.ref.csv',      'KELT-17b.median.csv'],$ ;; DT fitting
        [       'WASP-18b',  'fitwasp18',example_dir +      'wasp18' + path_sep(),       'WASP-18b.ref.csv',      'WASP-18b.median.csv'],$ ;; phase curve (+secondary) example
        [       'WASP-76b',  'fitwasp76',example_dir +      'wasp76' + path_sep(),       'WASP-76b.ref.csv',      'WASP-76b.median.csv'],$
;        [      'TOI-2165b', 'fittoi2165',example_dir +     'toi2165' + path_sep(),       'TOI2165b.ref.csv',     'TOI-2165b.median.csv'],$
        ['EPIC 247098361b',    'fit8361',example_dir + 'ep247098361' + path_sep(), 'EPIC247098361b.ref.csv','EPIC247098361b.median.csv'],$ ;; Kepler long cadence
        [         'GJ9827',  'fitgj9827',example_dir +      'gj9827' + path_sep(),         'GJ9827.ref.csv',        'GJ9827.median.csv'],$ ;; 3 planet, Kepler long cadence
;        [       'HIP50796',   'fit50796',example_dir +    'hip50796' + path_sep(),       'HIP50796.ref.csv',      'HIP50796.median.csv'],$ ;; astrometry...
        [       'HD106315','fithd106315',example_dir +    'hd106315' + path_sep(),       'HD106315.ref.csv',      'HD106315.median.csv']]  ;; multiplanet, transit only example
        
nfits = n_elements(fits[0,*])

for i=0L, nfits-1 do begin
   print, fits[0,i]
   if not keyword_set(compareonly) then begin
      call_procedure,fits[1,i], maxsteps=maxsteps, nthin=nthin, debug=debug, verbose=verbose, nthread=nthread
   endif

;   csvname = fits[2,i] + 'fitresults' + path_sep()+fits[4,i]
;   refname = fits[2,i] + fits[3,i]
;   comparecsv, csvname, refname, sigdiff=sigdiff, threshhold=threshhold
;   print
endfor

mkexofastv2

;; test license free uses (doesn't work on my local machine)
spawn, 'gdl -e "fithat3, maxsteps=' + strtrim(maxsteps,2) + '"'

argfile = filepath('hat3.args',$
                   root_dir=getenv('EXOFAST_PATH'),$
                   subdir=['examples','hat3_nolicense'])
vmfile = filepath('exofastv2.sav',$
                   root_dir=getenv('EXOFAST_PATH'))

spawn, 'idl -arg ' + argfile + ' -vm=' + vmfile

stop
end
