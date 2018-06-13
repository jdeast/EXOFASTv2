;+
; NAME:
;   exofast_errell
;
; PURPOSE: 
;   Gives the path of the error ellipses at constant probability given
;   two parameter distributions.
;
; PROCEDURE
;   1. Creates a 2D histogram of the data with NBINS bins per
;   parameter.
;   2. Calculates the fraction of data points above each value in the
;   histogram.
;   3. Compares each calculated fraction with the probabilities,
;   and defines the closest as the levels to pass to CONTOUR.
;   4. Creates a contour plot of the 2D histogram
;   5. Scales the path to be in the same units as the input and
;   recreates the path, adding NaNs between each level for easy
;   plotting
;
; CALLING SEQUENCE:
;   exofast_errell, xpar, ypar [,NXBIN=nxbin,NYBIN=nybin,PROBS=PROBS]$
;                   [,XPATH=xpath][,YPATH=ypath][,/PLOT][,/OPLOT]
;
; INPUTS:
;   XPAR       - A vector of parameters to be plotted on the x axis.
;   YPAR       - A vector of parameters to be plotted on the y axis.
;
; OPTIONAL INPUTS:
;   NXBIN      - The number of bins to take in each direction. Default
;                is 100. If you have only a small number of values for
;                xpar and ypar (10000), this should be smaller. 
;   NYBIN      - The number of bins to take in each direction. Default
;                is NXBIN.
;   XMIN/YMIN  - The minimum x/y value to make the contours
;   XMAX/YMAX  - The maximum x/y value to make the contours
;   PROBS      - An array containing the probability contours to
;                plot. The default is
;                probs = erf([1.d0,2.d0]/sqrt(2.d0)), or
;                probs = [0.6826894921370859d0,$
;                         0.9544997361036416d0]
;                Note: to have a meaningful contour enclosing 99.7% of
;                the points, you must have ~1 million points or make
;                NBINS smaller.
;  
; OPTIONAL KEYWORDS:
;   PLOT       - Plot the contours on a new plot
;   OPLOT      - Overplot the contours on an existing plot
;
; OPTIONAL OUTPUTS:
;   XPATH/YPATH - The x/y values of the contours. For more
;                 customizable plotting, return these and plot inside
;                 your own program. All contours can be plotted simply
;                 with IDL> plot, xpath, ypath
;                 NOTE: XPATH/YPATH contain NaNs between levels to make the
;                 plotting easier.
;   OUTSIDENDX  - The indices of the parameters outside of the
;                 outer-most contour (requires David Fanning's Inside
;                 procedure)
;
; EXAMPLE:
; plot the 68% and 95% contours, and every point outside them 
;
; x = randomn(seed,100000)
; y = randomn(seed,100000)
; exofast_errell, x,y,out=out,probs=[0.68,0.95],xpath=xpath,ypath=ypath
; plot, xpath, ypath, /iso
; if out[0] ne -1 then oplot, x[out],y[out],psym=3
;
; DEPENDENCIES:
; INSIDE (Coyote library):
; http://www.idlcoyote.com/programs/inside.pro
;
; NOTE: there is non-vectorized version of his code online, which will
; not work. Use the one at the link above.
;
; REVISION HISTORY:
;   2010/03/01 - Written: Jason Eastman - The Ohio State University
;   2010/11/05 - Added OUTSIDENDX output (inspired by Ben Shappee, OSU)
;                Added, NXBIN, XMIN, XMAX, NYBIN, YMIN, YMAX inputs.
;   2011/03/25 - Better memory management for INSIDE -- split up
;                array to check (takes a little longer, but actually
;                works with large arrays!)
;-
pro exofast_errell, xpar, ypar, probs=probs, xpath=xpath, ypath=ypath, $
                    plot=plot, oplot=oplot, outsidendx=outsidendx, $
                    nxbin=nxbin, xmin=xmin, xmax=xmax, $
                    nybin=nybin, ymin=ymin, ymax=ymax,logname=logname

;; 1 and 2 sigma contours (as defined in 1D)
if n_elements(probs) eq 0 then probs = erf([1.d0,2.d0]/sqrt(2.d0))

;; define the bounds of the histogram in the x direction
if n_elements(xmin) eq 0 then xmin = min(xpar,/nan)
if n_elements(xmax) eq 0 then xmax = max(xpar,/nan)
if n_elements(nxbin) eq 0 then nxbin = 100
if xmin eq xmax then begin
   printandlog, 'Histogram has zero width (parameter has no dispersion?)',logname
   return
endif
xbinsz = double(xmax-xmin)/(nxbin-1)

;; define the bounds of the histogram in the y direction
if n_elements(ymin) eq 0 then ymin = min(ypar,/nan)
if n_elements(ymax) eq 0 then ymax = max(ypar,/nan)
if n_elements(nybin) eq 0 then nybin = nxbin
if ymin eq ymax then begin
   printandlog, 'Histogram has zero width (parameter has no dispersion?)',logname
   return
endif
ybinsz = double(ymax-ymin)/(nybin-1)

;; create the 2D histogram
hist2d = hist_2d(xpar,ypar,$
                 bin1=xbinsz,min1=xmin,max1=xmax,$
                 bin2=ybinsz,min2=ymin,max2=ymax)


;; find the fraction of points enclosed above each bin 
;; NOTE: includes points excluded by xmin, xmax, ymin, ymax
maxhist = max(hist2d)
enclosed = dblarr(maxhist+1)
for i=0L, maxhist do begin
    good = where(hist2d ge i)
    if good[0] ne -1 then enclosed[i] = total(hist2d[good])
endfor
fraction = enclosed/n_elements(where(finite(xpar)))

;; find the levels that correspond to the different values of sigma
levels = dblarr(n_elements(probs))
for i=0, n_elements(probs)-1 do begin
    min = min(abs(fraction-probs[i]),level)
    levels[i] = level
endfor
levels = reverse(levels)

;; make sure LEVELS is a sensible input to CONTOUR
if n_elements(uniq(levels)) ne n_elements(probs) then begin
    printandlog, 'Resolution too fine -- lower NXBIN or use more elements!',logname
    return
endif

;; get the paths of the contours
contour, hist2d, xtitle=xtitle,ytitle=ytitle, levels=levels,$
  path_info=info,path_xy=xy,/path_data_coords,/overplot,/path_double

;; prep for outside index
if arg_present(outsidendx) then isinside = lonarr(n_elements(xpar))

;; convert the path values input units
x = xmin + (xy[0,*]+1)*(xmax-xmin)/nxbin
y = ymin + (xy[1,*]+1)*(ymax-ymin)/nybin

;; convert the path to intuitively plottable path
;; ie, plot, xpath, ypath
xpath = !values.d_nan
ypath = !values.d_nan
for i=0, n_elements(info)-1 do begin
    S = [indgen(info(i).n), 0]
    xpath = [xpath,x[info(i).offset + S],!values.d_nan]
    ypath = [ypath,y[info(i).offset + S],!values.d_nan]

    ;; finds elements that are outside all outer-level contours
    ;; (isinside = 0)
    if info(i).level eq 0 and arg_present(outsidendx) then begin
        check = where(isinside eq 0, ncheck)

        if check[0] ne -1 then begin        
            ;; inside is very memory inefficient --
            ;; split up the check to more managable bites
            for j=0L, ncheck, 1d4 do begin
                checknow = j + indgen(1d4 < (ncheck-j))
                isinside[check[checknow]] += $
                  inside(xpar[check[checknow]],ypar[check[checknow]],$
                         x[info(i).offset+S],y[info(i).offset+S])
            endfor 
        endif
    endif
endfor

if arg_present(outsidendx) then outsidendx = where(isinside eq 0)

;; begin a new plot if desired
if keyword_set(plot) then $
  plot, xpath, ypath,_extra=extra

;; overplot if desired
if keyword_set(oplot) then $
  oplot, xpath, ypath

end
