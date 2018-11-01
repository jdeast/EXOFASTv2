pro plotchain, pars, label=label, unit=unit, latex=latex, logname=logname, burnndx=burnndx

nallsteps = n_elements(pars[*,0])
nchains = n_elements(pars[0,*])
trimpars = pars[burnndx:*,*]
nsteps = n_elements(trimpars)
medndx = (nsteps/2d0)
halfsigma = erf(1d/sqrt(2d0))/2d
lowsigndx = round( nsteps/2.d0 - nsteps*halfsigma )
hisigndx = round( nsteps/2.d0 + nsteps*halfsigma )
charsize = 1.5

;; check for bad values
bad = where(~finite(pars),complement=good)
if bad[0] ne -1 then begin
   printandlog, "ERROR: NaNs in " + label + " distribution", logname
   return
endif

;; if angular, center distribution about the mode
if unit eq 'DEGREES' then halfrange = 180d0 $
else if unit eq 'RADIANS' then halfrange = !dpi
if unit eq 'DEGREES' or unit eq 'RADIANS' then begin
   
   ;; find the mode
   hist = histogram(pars,nbins=100,locations=x,/nan)
   max = max(hist,modendx)
   mode = x[modendx]
   
   toohigh = where(pars gt (mode + halfrange))
   if toohigh[0] ne -1 then pars[toohigh] -= 2.d0*halfrange
   
   toolow = where(pars lt (mode - halfrange))
   if toolow[0] ne -1 then pars[toolow] += 2.d0*halfrange
endif                     

;; 68% confidence interval
sorted = sort(trimpars)
medvalue = trimpars[sorted[medndx]]
upper = trimpars[sorted[hisigndx]] - medvalue
lower = medvalue - trimpars[sorted[lowsigndx]]

xmax = (medvalue + 4*upper) < max(pars)
xmin = (medvalue - 4*lower) > min(pars)

if xmin eq xmax then begin
   printandlog, 'WARNING: ' + label + ' is singularly valued.',logname
endif else begin
   
   ;; plot labels
   xtitle='!3 Chain link number'
   ytitle='!3' + exofast_textoidl(latex)
   
   plot, [0],[0], psym=10, xtitle=xtitle, ytitle=ytitle,$
         charsize=charsize,xstyle=1, yrange=[xmin,xmax], xrange=[0,nallsteps],font=1
   
   ;; thin the chain for plotting 
   ;; (otherwise the file becomes unusably large)
   nplot = 1d3 < n_elements(pars[*,0])
   thin = n_elements(pars[*,0])*lindgen(nplot)/(nplot)
   for l=0L, nchains-1L do $
      oplot, thin, pars[thin,l], color=l*255d0/nchains
   oplot, [burnndx,burnndx],[-9d99,9d99]
   
endelse

end

pro exofast_plotchains, ss, chainfile=chainfile, logname=logname

if n_elements(chainfile) eq 0 then chainfile = 'chain.ps'

chi2 = reform((*ss.chi2),ss.nsteps/ss.nchains,ss.nchains)
burnndx = getburnndx(chi2,goodchains=goodchains, logname=logname)

;; prepare the postscript device
mydevice=!d.name
set_plot, 'PS'
aspect_ratio=1.5
xsize = 18
ysize=xsize/aspect_ratio
!p.font=0

defsysv, '!GDL', exists=runninggdl  
if runninggdl then chainfile0 = file_dirname(chainfile) + path_sep() +  file_basename(chainfile,'.ps') + '.001.ps' $
else chainfile0 = chainfile

device, filename=chainfile0, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
loadct,39,/silent
red = 254

!p.multi=[0,2,4] ;; 8 to a page

ymin = min(chi2,max=ymax)
if ymin lt 0 then sign = ' + ' else sign = ' - '
latex = '\chi^2' + sign + strtrim(abs(ymin),2)
plotchain, chi2[*,goodchains]-ymin, latex=latex, unit='', label='chi2', logname=logname, burnndx=burnndx
page = 2

for i=0, n_tags(ss)-1 do begin
   for j=0, n_elements(ss.(i))-1 do begin
      for k=0, n_tags(ss.(i)[j])-1 do begin

         derive = 0B
         m=-1
         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin
            if ptr_valid(ss.(i)[j].(k)) then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'derive') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].derive or (*(ss.(i)[j].(k))).(l)[m].fit then begin
                              ;; GDL can't do multi-page plots
                              if runninggdl and !p.multi[0] eq 0 then begin
                                 device, /close
                                 chainfile0 = file_dirname(chainfile) + path_sep() + file_basename(chainfile,'.ps') + '.' + string(page,format='(i03)') + '.ps'
                                 device, filename=chainfile0, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
                                 page += 1
                              endif
                              pars = (reform((*(ss.(i)[j].(k))).(l)[m].value,ss.nsteps/ss.nchains,ss.nchains))[*,goodchains]
                              plotchain, pars, label=(*(ss.(i)[j].(k))).(l)[m].label,$
                                         unit = (*(ss.(i)[j].(k))).(l)[m].unit,$
                                         latex = (*(ss.(i)[j].(k))).(l)[m].latex,$
                                         logname=logname, burnndx=burnndx
                           endif
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            if n_tags(ss.(i)[j].(k)) ne 0 then begin
               if tag_exist(ss.(i)[j].(k),'derive') then begin
                  if ss.(i)[j].(k).derive then begin
                     ;; GDL can't do multi-page plots
                     if runninggdl and !p.multi[0] eq 0 then begin
                        device, /close
                        chainfile0 = file_dirname(chainfile) + path_sep() + file_basename(chainfile,'.ps') + '.' + string(page,format='(i03)') + '.ps'
                        device, filename=chainfile0, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
                        page += 1
                     endif
                     pars = (reform(ss.(i)[j].(k).value,ss.nsteps/ss.nchains,ss.nchains))[*,goodchains]
                     plotchain, pars, label=ss.(i)[j].(k).label,$ 
                                unit = ss.(i)[j].(k).unit,$
                                latex = ss.(i)[j].(k).latex,$
                                logname=logname, burnndx=burnndx
                  endif
               endif
            endif
         endif

      endfor
   endfor
endfor

;; clean up the postscript device
device, /close
!p.font=-1
!p.multi=0
set_plot, mydevice

end

