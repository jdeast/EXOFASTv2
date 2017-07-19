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
pro exofast_plotdist, pars, medianpars, angular=angular, $
                      parnames=parnames, units=units, nocovar=nocovar, $
                      pdfname=pdfname, covarname=covarname, $
                      bestpars=bestpars, probs=probs, degrees=degrees, logname=logname

if n_elements(pdfname) eq 0 then pdfname = 'pdf.ps'
if n_elements(covarname) eq 0 then covarname = 'covar.ps'
if n_elements(bestpars) ne 0 then bestpars = replicate(1,3)#bestpars

nparnames = n_elements(parnames)
nunits = n_elements(units)

;; 68% and 95% probability contours
if n_elements(probs) eq 0 then probs = erf([1d,2d]/sqrt(2d0))

sz = size(pars)
if sz[0] ne 2 then message, 'ERROR: PARS must be an NPARS x NSTEPS array'
npars = double(sz[1])
nsteps = double(sz[2])

if nparnames ne npars then begin
   parnames = strarr(npars)
   if nparnames ne 0 then printandlog, $
      'WARNING: PARNAMES does not match parameter array; ignoring labels',logname
endif

if nunits ne npars then begin
   units = strarr(npars)
   if nunits ne 0 then printandlog, $
      'WARNING: UNITS does not match parameter array; ignoring labels',logname
endif

medndx = nsteps/2
halfsigma = erf(1d/sqrt(2d0))/2d
lowsigndx = round(nsteps/2.d0 - nsteps*halfsigma)
hisigndx = round(nsteps/2.d0 + nsteps*halfsigma)
medianpars = dblarr(3,npars)

;; prepare the postscript device
mydevice=!d.name
set_plot, 'PS'
aspect_ratio=1.5
xsize = 18
ysize=xsize/aspect_ratio
!p.font=1
device, filename=pdfname
device, xsize=xsize,ysize=ysize

charsize=1.5

;; for angular parameters, center about the mode first
nangular = n_elements(angular)
for i=0, nangular-1 do begin

    hist = histogram(pars[angular[i],*],nbins=100,locations=x,/nan)
    max = max(hist,modendx)
    mode = x[modendx]

    ;; choose degrees vs radians
    if keyword_set(degrees) then halfrange = 180d0 $
    else halfrange = !dpi

    toohigh = where(pars[angular[i],*] gt (mode + halfrange))
    if toohigh[0] ne -1 then pars[angular[i],toohigh] -= 2.d0*halfrange

    toolow = where(pars[angular[i],*] lt (mode - halfrange))
    if toolow[0] ne -1 then pars[angular[i],toolow] += 2.d0*halfrange

endfor

!p.multi=[0,2,4] ;; 8 to a page
for i=0, npars-1 do begin

    ytitle='Probability'

    ;; create the x axis titles
    xtitle = parnames[i]
        
    ;; add the units
    if units[i] ne '' then xtitle += ' (' + units[i] + ')'
     
    bad = where(~finite(pars[i,*]),complement=good)
    if bad[0] ne -1 then printandlog, $
       "ERROR: NaNs in " +  parnames[i] + " distribution",logname

    ;; get the median and +/- 1 sigma values
    sorted = sort(pars[i,*])
    medianpars[0,i] = pars[i,sorted[medndx]]
    medianpars[1,i] = pars[i,sorted[hisigndx]] - medianpars[0,i] ; + error
    medianpars[2,i] = medianpars[0,i] - pars[i,sorted[lowsigndx]]; - error

    xmax = (medianpars[0,i] + 4*medianpars[1,i]) < max(pars[i,*],/nan)
    xmin = (medianpars[0,i] - 4*medianpars[2,i]) > min(pars[i,*],/nan)
    
    if xmin eq xmax then begin
       printandlog, 'WARNING: Parameter ' + strtrim(i,2) + ' (' + parnames[i] + $
          ') is singularly valued; this is often a symptom of a bigger problem. Skipping covariances, which would fail.',logname
       ;nocovar = 1
    endif

    ;; plot the probability distribution
    hist = histogram(pars[i,sorted],nbins=100,locations=x,min=xmin,max=xmax,/nan)
    plot, x, hist/double(nsteps), psym=10, xtitle=xtitle, ytitle=ytitle,$
      charsize=charsize,xstyle=1, xrange=[xmin,xmax]
    
    ;; if the best parameters are given, overplot them on the PDFs and
    ;; give the confidence interval around them
    if n_elements(bestpars) eq 3*npars then begin
       oplot, [bestpars[0,i],bestpars[0,i]],[-9d9,9d9]

       dummy = min(abs(pars[i,sorted]-bestpars[0,i]),bestndx) 
       lowndx = round(bestndx - nsteps*halfsigma)
       hindx = round(bestndx + nsteps*halfsigma)

       ;; + error (for skewed distributions, it's possible this
       ;; doesn't exist; make it infinity)
       if (hindx le nsteps-1) then $
          bestpars[1,i] = pars[i,sorted[hindx]] - bestpars[0,i] $
       else bestpars[1,i] = !values.d_infinity

       ;; - error (for skewed distributions, it's possible this
       ;; doesn't exist; make it infinity)
       if lowndx ge 0 then $
          bestpars[2,i]=bestpars[0,i] - pars[i,sorted[lowndx]] $
       else bestpars[2,i] = !values.d_infinity
       
    endif

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
nplot=0
nplots = npars*(npars-1.d0)/2.d0
nxplots = 4                     
nyplots = 4   
corr = correlate(pars)

!p.multi=[0,nxplots,nyplots]
for i=0, npars-1 do begin
   xtitle = parnames[i]
   xmax = (medianpars[0,i] + 4*medianpars[1,i]) < max(pars[i,*],/nan)
   xmin = (medianpars[0,i] - 4*medianpars[2,i]) > min(pars[i,*],/nan)

   for j=i+1,npars-1 do begin
      
      ;; create the titles
      ytitle = parnames[j]
      
      ;; units just clutter the page, not really necessary
      ;if units[i] ne '' then xtitle += ' (' + units[i] + ')'
      ;if units[j] ne '' then ytitle += ' (' + units[j] + ')'

      ;; the ranges of the plots/histograms
      ymax = (medianpars[0,j] + 4*medianpars[1,j]) < max(pars[j,*],/nan)
      ymin = (medianpars[0,j] - 4*medianpars[2,j]) > min(pars[j,*],/nan)
      
      if xmin ne xmax and ymin ne ymax then begin

         ;; get the paths of the contours
         exofast_errell, pars[i,*],pars[j,*],xpath=xpath,ypath=ypath,$
                         prob=probs,nxbin=100, nybin=100,$
                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax
         
         ;; plot the contours
         plot, xpath, ypath, xtitle=xtitle,ytitle=ytitle, $
               xrange=[xmin,xmax],yrange=[ymin,ymax],/ystyle,/xstyle,$
               xtickname=replicate(" ",60),ytickname=replicate(" ",60),$
               charsize=1.5,title=strtrim(corr[i,j],2)
         
         ;; plot the best fit lines
         if n_elements(bestpars) eq npars*3 then begin
            plotsym, 0, 0.5, /fill 
            oplot, [bestpars[0,i]],[bestpars[0,j]], psym=8
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
