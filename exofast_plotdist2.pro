;+
; NAME:
;   exofast_plotdist
; PURPOSE:
;   Plots the distributions from MCMC chains (Probability Distribution
;   Functions (PDFs) and covariances; returns median values and 68%
;   confidence intervals.
;
; DESCRIPTION:
;
; CALLING SEQUENCE:
;   exofast_plotdist, pars [,MEDIANPARS, ANGULAR=, PARNAMES=, UNITS=,
;                     /NOCOVAR, PDFNAME=, COVARNAME=, BESTPARS=, PROBS=, /DEGREES]
;
; INPUTS:
;   PARS       - Array of NPARS x NSTEPS for each step in a Markov
;                Chain (with all chains concatenated together).
;
; OPTIONAL INPUTS:
;   ANGULAR    - An array that indexes the parameters that are
;                angular (radians unless /DEGREES is set). This will enable
;                special ways to calculate the median and plots, which
;                may otherwise fail. Default is none.
;   PARNAMES   - An NPARS string array containing the parameter names
;                to be plotted on the axes.
;                HINT: TeXtoIDL is awesome.
;   UNITS      - An NPARS string array containing the units to be
;                displayed on the X-axis of the PDFs.
;   PDFNAME    - The filename of the output postscript plot of
;                PDFs. Default is pdf.ps
;   COVARNAME  - The filename of the output postscript plot of
;                covariances. Default is covar.ps. 
;   BESTPARS   - An NPARS array that specifies the "best" parameters
;                of each value. If specified, a line will be drawn on
;                each PDF and a dot will be drawn on each covariance
;                plot at these values (this doesn't have to be
;                the best values).
;   PROBS      - Array of probability contours to draw on the covariance
;                plots. Default is erf([1d,2d]/sqrt(2d0)) ~ [0.68,0.95].
;
; OPTIONAL KEYWORDS:
;   DEGREES    - If set, ANGULAR parameters are in degrees, not radians.
;   NOCOVAR    - If set, no covariances plots will be made. Since
;                there will be NPARS*(NPARS-1)/2 plots, this can be
;                time consuming.
; OUTPUTS:
;   MEDIANPARS - A 3 x NPARS array of median values, +err, -err (68%
;                confidence interval)
;   PDFNAME    - A postscript plot of Probablity Distribution
;                Functions (PDFs)
;   COVARNAME  - A postscript plot of covariances.
;
; REVISION HISTORY:
;   2009/11/24 - Written: Jason Eastman - The Ohio State University
;   2013/03/12 - Fixed min/max, made slightly more efficient
;-
pro exofast_plotdist2, ss, nocovar=nocovar, $
                       pdfname=pdfname, covarname=covarname, $
                       probs=probs,logname=logname

if n_elements(pdfname) eq 0 then pdfname = 'pdf.ps'
if n_elements(covarname) eq 0 then covarname = 'covar.ps'

;; 68% and 95% probability contours
if n_elements(probs) eq 0 then probs = erf([1d,2d]/sqrt(2d0))

burnndx = ss.burnndx
nsteps = ss.nsteps/ss.nchains
medndx = ((nsteps-ss.burnndx)/2d0) * ss.nchains
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
device, filename=pdfname
device, set_font='Times',/tt_font
device, xsize=xsize,ysize=ysize,/color,bits=24
charsize=1.5
loadct,39,/silent
red = 254

allpars = []
parnames = []
medianpars = []

!p.multi=[0,2,4] ;; 8 to a page
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
                              pars = (reform((*(ss.(i)[j].(k))).(l)[m].value,nsteps,ss.nchains))[burnndx:*,*]
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
                     pars = (reform(ss.(i)[j].(k).value,nsteps,ss.nchains))[burnndx:*,*]
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
            if bad[0] ne -1 then printandlog, $
               "ERROR: NaNs in " + label + " distribution",logname
            
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
            xmin = (medvalue - 4*lower) > min(pars)
            parnames = [parnames,latex]
            
            if xmin eq xmax then begin
               printandlog, 'WARNING: ' + label + ' is singularly valued.',logname
            endif else begin
               
               ;; plot labels
               xtitle='!3' + textoidl(latex + '_{' + ss.(i)[j].label + '}')
               ytitle='!3Probability'
              


               hist = histogram(pars,nbins=100,locations=x,min=xmin,max=xmax)
               plot, x, hist/double(total(hist)), psym=10, xtitle=xtitle, ytitle=ytitle,$
                     charsize=charsize,xstyle=1, xrange=[xmin,xmax], font=1

               for l=0L, ss.nchains-1L do  begin
                  hist = histogram(pars[*,l],nbins=100,locations=x,min=xmin,max=xmax)
                  if total(hist) gt 0 then $
                     oplot, x, hist/double(total(hist)), psym=10,color=l*255d0/ss.nchains;, transparency=100d0-100d0/ss.nchains

               endfor

               hist = histogram(pars,nbins=100,locations=x,min=xmin,max=xmax)
               oplot, x, hist/double(total(hist)), psym=10, thick=3
               
               ;; if the best parameters are given, overplot them on the PDFs
               if finite(best) then $
                  oplot, [best,best],[-9d9,9d9]
            endelse
            
            
            ;; format values for table (rounded appropriately)
            ;; round the high error to 2 sig figs
            exphi=fix(alog10(upper))
            if (upper lt 1d0) then exphi=exphi-1
            roundhi=round(upper/10.d0^(exphi-1d0),/L64)*10.d0^(exphi-1d0)
            if (roundhi gt 10) then errhi = strtrim(round(roundhi,/L64),2) $
            else errhi = string(roundhi,format='(f255.'+strtrim(1-exphi,2)+')')
            
            ;; round the low error to 2 sig figs
            explo=fix(alog10(lower))
            if (lower lt 1d0) then explo=explo-1
            roundlo=round(lower/10.d0^(explo-1d0),/L64)*10.d0^(explo-1d0)
            if (roundlo gt 10) then errlo = strtrim(round(roundlo,/L64),2) $
            else errlo = string(roundlo,format='(f255.'+strtrim(1-explo,2)+')')
            
            ;; round the value to the greater number of sig figs
            ndec = long(1 - (exphi < explo))
            if ndec eq 0 then value = string(medvalue,format='(i255)') $
            else if ndec lt 0 then $
               value=round(round(medvalue/10.d0^(-ndec),/L64)*10.d0^(-ndec),/L64) $
            else value = string(medvalue,format='(f255.'+strtrim(ndec,2)+')')
            
            if m ne -1 then begin
               (*ss.(i)[j].(k)).(l)[m].medvalue = strtrim(value,2)
               (*ss.(i)[j].(k)).(l)[m].upper = strtrim(errhi,2)
               (*ss.(i)[j].(k)).(l)[m].lower = strtrim(errlo,2)
            endif else begin
               ss.(i)[j].(k).medvalue = strtrim(value,2)
               ss.(i)[j].(k).upper = strtrim(errhi,2)
               ss.(i)[j].(k).lower = strtrim(errlo,2)
            endelse

         endif
      endfor
   endfor
endfor

device, /close

;; if covariance plots aren't wanted (they take a while), we're done
if keyword_set(nocovar) then begin
   !p.font=-1
   !p.multi=0
   set_plot, mydevice
   return
endif

device, filename = covarname
nxplots = 4                     
nyplots = 4   
corr = correlate(allpars)
npars = n_elements(allpars[*,0])
nplot = 0

!p.multi=[0,nxplots,nyplots]
for i=0, npars-1 do begin
   xtitle = '!3'+parnames[i]
   xmax = (medianpars[0,i] + 4*medianpars[1,i]) < max(allpars[i,*],/nan)
   xmin = (medianpars[0,i] - 4*medianpars[2,i]) > min(allpars[i,*],/nan)

   for j=i+1,npars-1 do begin
      
      ;; create the titles
      ytitle = '!3'+parnames[j]
      
      ;; units just clutter the page, not really necessary
      ;if units[i] ne '' then xtitle += ' (' + units[i] + ')'
      ;if units[j] ne '' then ytitle += ' (' + units[j] + ')'

      ;; the ranges of the plots/histograms
      ymax = (medianpars[0,j] + 4*medianpars[1,j]) < max(allpars[j,*],/nan)
      ymin = (medianpars[0,j] - 4*medianpars[2,j]) > min(allpars[j,*],/nan)
      
      if xmin ne xmax and ymin ne ymax then begin

         xpath = []
         ypath = []

         ;; get the paths of the contours
         exofast_errell, allpars[i,*],allpars[j,*],xpath=xpath,ypath=ypath,$
                         prob=probs,nxbin=100, nybin=100,$
                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,logname=logname
         
         if n_elements(xpath) gt 0 and n_elements(ypath) eq n_elements(xpath) then begin
            ;; plot the contours
            plot, xpath, ypath, xtitle=xtitle,ytitle=ytitle, $
                  xrange=[xmin,xmax],yrange=[ymin,ymax],/ystyle,/xstyle,$
                  xtickname=replicate(" ",60),ytickname=replicate(" ",60),$
                  charsize=1.5,title=strtrim(corr[i,j],2),font=1
         
            ;; plot the best fit lines
            if n_elements(bestpars) eq npars*3 then begin
               plotsym, 0, 0.5, /fill 
               oplot, [bestpars[0,i]],[bestpars[0,j]], psym=8
            endif
         endif
      endif
      nplot++
      
   endfor
endfor

;; clean up the postscript device
device, /close
!p.font=-1
!p.multi=0
set_plot, mydevice

end
