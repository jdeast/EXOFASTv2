pro exofast_plotchains, ss, chainfile=chainfile

if n_elements(chainfile) eq 0 then chainfile = 'chain.ps'

burnndx = ss.burnndx
nsteps = ss.nsteps/ss.nchains
medndx = (nsteps-ss.burnndx)/2
halfsigma = erf(1d/sqrt(2d0))/2d
lowsigndx = round(((nsteps-ss.burnndx)/2.d0 - (nsteps-ss.burnndx)*halfsigma)*ss.nchains)
hisigndx = round(((nsteps-ss.burnndx)/2.d0 + (nsteps-ss.burnndx)*halfsigma)*ss.nchains)
;; prepare the postscript device
mydevice=!d.name
set_plot, 'PS'
aspect_ratio=1.5
xsize = 18
ysize=xsize/aspect_ratio
!p.font=0
device, filename=chainfile
device, set_font='Times',/tt_font
device, xsize=xsize,ysize=ysize,/color,bits=24
charsize=1.5
loadct,39,/silent
red = 254

allpars = []
parnames = []
medianpars = []
!p.multi=[0,2,4] ;; 8 to a page

chi2 = reform((*ss.chi2),ss.nsteps/ss.nchains,ss.nchains)
ymin = min(chi2[ss.burnndx:-1,*],max=ymax)
xtitle='!3 Chain link number'
if ymin lt 0 then sign = ' + ' else sign = ' - '
ytitle='!3' + textoidl('\chi^2') + sign + strtrim(abs(ymin),2)
plot, [0],[0], psym=10, xtitle=xtitle, ytitle=ytitle,$
      charsize=charsize,xstyle=1, yrange=[0,ymax-ymin], xrange=[0,nsteps],font=1
for l=0L, ss.nchains-1L do $
   oplot, chi2[*,l] - ymin, color=l*255d0/ss.nchains ;, transparency=100d0-100d0/ss.nchains
oplot, [ss.burnndx,ss.burnndx],[-9d99,9d99]

for i=0, n_tags(ss)-1 do begin
   for j=0, n_elements(ss.(i))-1 do begin
      for k=0, n_tags(ss.(i)[j])-1 do begin

         derive = 0B
         m=-1
         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin
            if (size(ss.(i)[j].(k)))[0] ne 0 then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'derive') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].derive then begin
                              pars = (reform((*(ss.(i)[j].(k))).(l)[m].value,nsteps,ss.nchains))
;                              pars = (*(ss.(i)[j].(k))).(l)[m].value[burnndx:*,*]
                              derive = 1B
                              label = (*(ss.(i)[j].(k))).(l)[m].label
                              unit = (*(ss.(i)[j].(k))).(l)[m].unit
                              latex = (*(ss.(i)[j].(k))).(l)[m].latex
                              best = (*(ss.(i)[j].(k))).(l)[m].best
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
                     pars = (reform(ss.(i)[j].(k).value,nsteps,ss.nchains))
                     derive = 1B
                     label = ss.(i)[j].(k).label 
                     unit = ss.(i)[j].(k).unit
                     latex = ss.(i)[j].(k).latex
                     best = ss.(i)[j].(k).best
                  endif
               endif
            endif
         endif

         if derive then begin

            if n_elements(allpars) eq 0 then allpars = transpose(pars) $
            else allpars = [allpars,transpose(pars)]

            ;; check for bad values
            bad = where(~finite(pars),complement=good)
            if bad[0] ne -1 then message, $
               "ERROR: NaNs in " + label + " distribution"
            
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
            sorted = sort(pars)
            medvalue = pars[sorted[medndx]]
            upper = pars[sorted[hisigndx]] - medvalue
            lower = medvalue - pars[sorted[lowsigndx]]
            
            if n_elements(medianpars) eq 0 then medianpars = [medvalue,upper,lower] $
            else medianpars = [[medianpars],[[medvalue,upper,lower]]]
            
            xmax = (medvalue + 4*upper) < max(pars)
            xmin = (medvalue + 4*lower) > min(pars)
            parnames = [parnames,latex]
            
            if xmin eq xmax then begin
               message, 'WARNING: ' + label + ' is singularly valued.',/continue
            endif else begin
               
               ;; plot labels
               xtitle='!3 Chain link number'
               ytitle='!3' + textoidl(latex + '_{' + ss.(i)[j].label + '}')
                             
               plot, [0],[0], psym=10, xtitle=xtitle, ytitle=ytitle,$
                     charsize=charsize,xstyle=1, yrange=[xmin,xmax], xrange=[0,nsteps],font=1

               for l=0L, ss.nchains-1L do $
                  oplot, pars[*,l], color=l*255d0/ss.nchains ;, transparency=100d0-100d0/ss.nchains
               oplot, [ss.burnndx,ss.burnndx],[-9d99,9d99]

            endelse
         endif
      endfor
   endfor
endfor

end
