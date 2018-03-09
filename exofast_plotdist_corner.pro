;+
; NAME:
;   exofast_plotdist_corner
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
;   SS         - The stellar structure created by EXOFASTv2, populated
;                with an array of values for each parameter
;
; OPTIONAL INPUTS:
;   ANGULAR    - An array that indexes the parameters that are
;                angular (radians unless /DEGREES is set). This will enable
;                special ways to calculate the median and plots, which
;                may otherwise fail. Default is none.
;   PDFNAME    - The filename of the output postscript plot of
;                PDFs. Default is pdf.ps
;   COVARNAME  - The filename of the output postscript plot of
;                covariances of fitted parameters. Default is covar.ps. 
;   PROBS      - Array of probability contours to draw on the covariance
;                plots. Default is erf([1d,2d]/sqrt(2d0)) ~ [0.68,0.95].
;   LOGNAME    - The name of the logfile to print messages to.
;
; OPTIONAL KEYWORDS:
;   NOCOVAR    - If set, no covariances plots will be made. Since
;                there will be NFIT*(NFIT-1)/2 plots, this can be
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
;   2018/02/13 - Covariance plot now a corner plot inspired by
;                corner.py (much prettier)
;-
pro exofast_plotdist_corner, ss, nocovar=nocovar, $
                             pdfname=pdfname, covarname=covarname, $
                             probs=probs,logname=logname, angular=angular,$
                             useallpars=useallpars, csvfile=csvfile

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

allpars = [] ;; a 2D array of all parameters for the covariance matrix
allparndx = [] ;; an map from the median values of all parameters to the covariance parameters
bestpars = [] ;; the best-fit (most likely) values to overplot on the PDF and covariance plots
parnames = [] ;; parameter names for the axis labels
medianpars = [] ;; median parameters and 68% confidence intervals
alltitles=[]

minchi2 = min(*ss.chi2,bestndx)

!p.multi=[0,2,4] ;; 8 to a page
for i=0, n_tags(ss)-1 do begin
   for j=0, n_elements(ss.(i))-1 do begin
      for k=0, n_tags(ss.(i)[j])-1 do begin

         derive = 0B
         fit = 0B
         m=-1
         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin
            if (size(ss.(i)[j].(k)))[0] ne 0 then begin
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin 
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'derive') then begin
                           if (*(ss.(i)[j].(k))).(l)[m].derive or (*(ss.(i)[j].(k))).(l)[m].fit then begin
                              pars = (reform((*(ss.(i)[j].(k))).(l)[m].value,nsteps,ss.nchains))[burnndx:*,*]
                              derive = 1B
                              if (*(ss.(i)[j].(k))).(l)[m].fit then fit = 1B
                              label = (*(ss.(i)[j].(k))).(l)[m].label
                              unit = (*(ss.(i)[j].(k))).(l)[m].unit
                              latex = (*(ss.(i)[j].(k))).(l)[m].latex
                              best = (*(ss.(i)[j].(k))).(l)[m].best
                              if ~finite(best) then best = (*(ss.(i)[j].(k))).(l)[m].value[bestndx]
                           endif
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            if n_tags(ss.(i)[j].(k)) ne 0 then begin
               if tag_exist(ss.(i)[j].(k),'derive') then begin
                  if ss.(i)[j].(k).derive or ss.(i)[j].(k).fit then begin
                     pars = (reform(ss.(i)[j].(k).value,nsteps,ss.nchains))[burnndx:*,*]
                     derive = 1B                     
                     if ss.(i)[j].(k).fit then fit = 1B
                     label = ss.(i)[j].(k).label 
                     unit = ss.(i)[j].(k).unit
                     latex = ss.(i)[j].(k).latex
                     best = ss.(i)[j].(k).best
                     if ~finite(best) then best = ss.(i)[j].(k).value[bestndx]
                  endif
               endif
            endif
         endif

         if fit or (keyword_set(useallpars) and derive) then begin
            sz = size(pars)
            if n_elements(allpars) eq 0 then allpars = transpose(reform(pars,sz[1]*sz[2])) $
            else allpars = [allpars,transpose(reform(pars,sz[1]*sz[2]))]
            parnames = [parnames,latex]
            if n_elements(medianpars) eq 0 then allparndx = [allparndx,0] $
            else allparndx = [allparndx,n_elements(medianpars[0,*])]
            bestpars = [bestpars,best]
         endif

         if derive then begin

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
            
            if xmin eq xmax then begin
               printandlog, 'WARNING: ' + label + ' is singularly valued.',logname
            endif else begin
               
               ;; plot labels
               xtitle='!3' + exofast_textoidl(latex + '_{' + ss.(i)[j].label + '}')
               ytitle='!3Probability'
              
               hist = histogram(pars,nbins=100,locations=x,min=xmin,max=xmax)
               plot, x, hist/double(total(hist)), psym=10, xtitle=xtitle, ytitle=ytitle,$
                     charsize=charsize,xstyle=1, xrange=[xmin,xmax], font=1

               for l=0L, ss.nchains-1L do  begin
                  hist = histogram(pars[*,l],nbins=100,locations=x,min=xmin,max=xmax)
                  if total(hist) gt 0 then $
                     oplot, x, hist/double(total(hist)), psym=10,color=l*255d0/ss.nchains
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
            
            alltitles = [alltitles,strtrim(value,2) + '^{+' + strtrim(errhi,2) + '}_{-' + strtrim(errlo,2) + '}']
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
;corr = correlate(allpars)
npars = n_elements(allpars[*,0])

!p.multi=0
charsize = 4d0/npars
multiplot, [npars,npars], /square, mXtitle='', mYtitle='', /rowmajor
for i=0, npars-1 do begin
  
   ;; x range of the plots/histograms
   xmax = (medianpars[0,allparndx[i]] + 4*medianpars[1,allparndx[i]]) < max(allpars[i,*],/nan)
   xmin = (medianpars[0,allparndx[i]] - 4*medianpars[2,allparndx[i]]) > min(allpars[i,*],/nan)

   for j=0,npars-1 do begin

      ;; the axis labels
      if i eq 0 then ytitle='!3' + exofast_textoidl(parnames[j]) else ytitle = ''
      if j eq npars-1 then xtitle='!3' + exofast_textoidl(parnames[i]) else xtitle = ''

      ;; skip the other half of the corner
      if j lt i+1 then begin

         ;; make a histogram in the diagonal elements
         if j eq i then begin
            hist = histogram(allpars[i,*],locations=x,nbin=100,min=xmin,max=xmax)
            plot, x, hist/total(hist), psym=10,xminor=1, yminor=1,$ 
                  xtickname=replicate(" ",60),ytickname=replicate(" ",60),xtitle=xtitle,$
                  charsize=charsize,xrange=[xmin,xmax],/xstyle
            
            ;; plot lines denoting the median and 68% confidence interval
            median = medianpars[0,allparndx[i]]
            losig = medianpars[0,allparndx[i]] - medianpars[1,allparndx[i]]
            hisig = medianpars[0,allparndx[i]] + medianpars[1,allparndx[i]]
            oplot, [median,median],[-9d9,9d9]
            oplot, [losig,losig],[-9d9,9d9],linestyle=1
            oplot, [hisig,hisig],[-9d9,9d9],linestyle=1

         endif

         multiplot ;; next plot
         continue
      endif
      
      ;; y range of the plots/histograms
      ymax = (medianpars[0,allparndx[j]] + 4*medianpars[1,allparndx[j]]) < max(allpars[j,*],/nan)
      ymin = (medianpars[0,allparndx[j]] - 4*medianpars[2,allparndx[j]]) > min(allpars[j,*],/nan)
      
      if xmin ne xmax and ymin ne ymax then begin

         xpath = []
         ypath = []

         ;; get the paths of the 68% and 95% contours
         exofast_errell, transpose(allpars[i,*]),transpose(allpars[j,*]),xpath=xpath,ypath=ypath,$
                         prob=probs,nxbin=100, nybin=100,$
                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,logname=logname
         
         ;; plot the contours if it worked
         if n_elements(xpath) gt 0 and n_elements(ypath) eq n_elements(xpath) then begin

            plot, xpath, ypath, xtitle=xtitle,ytitle=ytitle,charsize=charsize, $
                  xrange=[xmin,xmax],yrange=[ymin,ymax],/ystyle,/xstyle,xminor=1, yminor=1,$ 
                  xtickname=replicate(" ",60),ytickname=replicate(" ",60)
         
            ;; overplot the best fit point
            if n_elements(bestpars) eq npars then begin
               plotsym, 0, 3.5d0/npars, /fill
               oplot, [medianpars[0,allparndx[i]]],[medianpars[0,allparndx[j]]], psym=8
            endif

         endif
      endif
      multiplot ;; next plot
   endfor

endfor

;; reset back to default plotting parameters
multiplot, /reset
multiplot, /default 

;; clean up the postscript device
device, /close
!p.font=-1
!p.multi=0
set_plot, mydevice

end
