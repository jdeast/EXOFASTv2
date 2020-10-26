pro getcovarmatrix, filename, parnames=parnames, corr=corr, useallpars=useallpars, redo=redo

dirname = file_dirname(filename)
outname = filepath(file_basename(filename,'.mcmc.idl') + '.covar.idl',root_dir=dirname)

if file_test(outname) and ~keyword_set(redo) then return

restore, filename

if n_elements(pdfname) eq 0 then pdfname = 'pdf.ps'
if n_elements(covarname) eq 0 then covarname = 'covar.ps'

;; 68% and 95% probability contours
if n_elements(probs) eq 0 then probs = erf([1d,2d]/sqrt(2d0))

nsteps = mcmcss.nsteps/mcmcss.nchains
chi2 = reform((*mcmcss.chi2),nsteps,mcmcss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains)
ngoodchains = n_elements(goodchains)

allpars = [] ;; a 2D array of all parameters for the covariance matrix
parnames = [] ;; parameter names for the axis labels

for i=0, n_tags(mcmcss)-1 do begin
   for j=0, n_elements(mcmcss.(i))-1 do begin
      for k=0, n_tags(mcmcss.(i)[j])-1 do begin

         ;; this captures the detrending variables
         if (size(mcmcss.(i)[j].(k)))[1] eq 10 then begin ;; if it's a pointer
            if ptr_valid(mcmcss.(i)[j].(k)) then begin ;; if it's not empty
               for l=0L, n_tags(*(mcmcss.(i)[j].(k)))-1 do begin ;; loop through each tag
                  if (size((*(mcmcss.(i)[j].(k))).(l)))[2] eq 8 then begin ;; if it's an array of structures
                     for m=0L, n_elements((*(mcmcss.(i)[j].(k))).(l))-1 do begin ;; loop through each structure
                        if tag_exist((*(mcmcss.(i)[j].(k))).(l)[m],'derive') then begin ;; if it's a parameter
                           if (*(mcmcss.(i)[j].(k))).(l)[m].derive or (*(mcmcss.(i)[j].(k))).(l)[m].fit then begin ;; if it's requested to be fit or 
                              
                              ;; remove the burn-in, discard bad chains
                              pars = (reform((*(mcmcss.(i)[j].(k))).(l)[m].value,nsteps,mcmcss.nchains))[burnndx:*,goodchains]

                              ;; store these for the covariance plot
                              if (*(mcmcss.(i)[j].(k))).(l)[m].fit or keyword_set(useallpars) then begin
                                 sz = size(pars)
                                 if n_elements(allpars) eq 0 then allpars = transpose(reform(pars,sz[1]*sz[2])) $
                                 else allpars = [allpars,transpose(reform(pars,sz[1]*sz[2]))]
                                 parnames = [parnames,(*(mcmcss.(i)[j].(k))).(l)[m].latex]
                              endif
                     
                           endif
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(mcmcss.(i)[j].(k)) ne 0 then begin
            ;; and all other parameters
            if n_tags(mcmcss.(i)[j].(k)) ne 0 then begin
               if tag_exist(mcmcss.(i)[j].(k),'derive') then begin
                  if mcmcss.(i)[j].(k).derive or mcmcss.(i)[j].(k).fit then begin

                     ;; remove the burn-in, discard bad chains
                     pars = (reform(mcmcss.(i)[j].(k).value,nsteps,mcmcss.nchains))[burnndx:*,goodchains]

                     ;; store these for the covariance plot
                     if mcmcss.(i)[j].(k).fit or keyword_set(useallpars) then begin
                        sz = size(pars)
                        if n_elements(allpars) eq 0 then allpars = transpose(reform(pars,sz[1]*sz[2])) $
                        else allpars = [allpars,transpose(reform(pars,sz[1]*sz[2]))]
                        parnames = [parnames,mcmcss.(i)[j].(k).latex]
                     endif

                  endif
               endif
            endif
         endif
      endfor
   endfor
endfor

corr = correlate(allpars)
save, parnames, corr, filename=outname

end
