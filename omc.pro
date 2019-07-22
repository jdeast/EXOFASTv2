;+
; NAME:
;   OMC
;
; PURPOSE: 
;   Calculates the best ephemeris from a list of transit times and
;   plots the O-C diagram
; 
; CALLING SEQUENCE:
;    OMC, 'filename' [,period, t0=t0,/ps,epsname=epsname]
;
; INPUTS:
;
;    FILENAME - The filename of a text file with columns T_C, error,
;               and optionally telescope name. If telescope name is
;               given, each unique telescope name will be plotted
;               with a different color. The error must be in days.
;               NOTE: "#" may be used as a comment string.
;
; OPTIONAL INPUTS:
;    
;    PERIOD   - The period of the transit. If not given, the program
;               will attempt to calculate it from the lowest common
;               multiple of transits.
;               NOTE: This may result in integer aliases depending on
;               when the transits were taken.
;    T0       - The zero point of the ephemeris. If not specified, the
;               transit time closest to the error-weighted mean of all
;               input times will be used, which minimizes the error in
;               the ephemeris.
;    EPSNAME  - The name of the encapsulated postscript file; default
;               is omc.eps.
;
; OPTIONAL KEYWORDS:
;    
;    PS       - If set, an encapsulated postscript plot will be
;               created. Otherwise, it will plot to the screen. It is
;               not necessary to specify both PS and EPSNAME.
;
; OUTPUTS:
;
;   1) Plots the O-C diagram
;   2) Prints the best fit ephemeris with errors
;   3) Prints the chi^2/dof of the best fit
;   4) Prints the epoch, transit time, error (sec), residual (sec),
;      and residual/err for each transit
;
; MODIFICATION HISTORY
;
;  2011/07 -- Written, Jason Eastman (OSU)
;
;-
pro omc, time, err, telescope=telescope, chi2=chi2, period=period, t0=t0, ps=ps, epsname=epsname

;; read in the times
;openr, lun, filename, /get_lun
;line = ""
;readf, lun, line
;ncols = n_elements(strsplit(line))
;free_lun, lun
;if ncols eq 2 then begin
;    readcol, filename, time, err, format='d,d', comment='#',/silent 
;endif else if ncols eq 3 then begin
;    readcol, filename, time,err,telescope, format='d,d,a', comment='#',/silent
;endif else message, 'FILENAME must have 2 or three columns'

;; Determine the lowest common multiple of transit times 
;; this is either the period or an integer multiple of it
if n_elements(period) eq 0 then begin
    niter = 1d0
    diff = (time-shift(time,1))
    good = where(diff gt 0.1) ;; if diff lt 0.1, it's the same transit
    mindiff = min(diff[good])
    repeat begin        
        trialP = mindiff/niter
        np = diff[good]/trialP
        niter++
        if niter ge 10 then $
          message, 'Cannot determine best period; please input manually'
    endrep until max(abs(np - round(np))) lt 0.1
    period = trialP
endif

;; find the best t0 (i.e., nearest to the error-weighted mean of T_Cs)
if n_elements(t0) eq 0 then begin
    tmean = total(time*err)/total(err)
    np = round((time[0] - tmean)/period)
    t0 = time[0] - np*period
endif

;; calculate the epoch numbers (x axis)
epoch = round((time - t0)/period)

;; fit a straight line to the Transit Times
;; TC = coeffs[0] +/- sigma[0]
;; P  = coeffs[1] +/- sigma[1]
coeffs = poly_fit(epoch,time,1,measure_errors=err,yfit=yfit,sigma=sigma)

;; calculate the residuals
omc = (time-yfit)*1440d0
err *= 1440d0 ;; convert error to minutes

;; prepare the plotting device (thanks astrobetter)
;; http://www.astrobetter.com/making-fonts-better-in-idl-postscript-output/
if keyword_set(ps) or n_elements(epsname) ne 0 then begin
    if n_elements(epsname) eq 0 then epsname = 'omc.eps'
    mydevice = !d.name
    set_plot, 'PS'
    !p.font=0
    aspect_ratio=1.5
    xsize=9
    ysize=xsize/aspect_ratio
    device,  filename=epsname,/color,bits=24,encapsulated=1, /helvetica
    device, xsize=xsize, ysize=ysize
    LOADCT, 39,/silent
    colors = [0,254,159,95,223,31,207,111,191,47]
    charsizelegend = 0.09
    ;; left
    xlegend = 0.1

    ;; top
    ylegend = 0.90

    ;; right
;    xlegend = 0.80

    ;; bottom
;    ylegend = 0.60


    charsize = 0.5
endif else begin
    colors= ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x,'0080ff'x,$
             '800080'x,'00ffff'x,'ffff00'x,'80d000'x,'660000'x]
    charsizelegend = 0.03
    xlegend = 0.90
    ylegend = 0.99
    charsize = 1
endelse
ncolors = n_elements(colors)

;; plot the times, errors, and 0 line
ymin = min(omc-err) & ymax = max(omc+err)
xmin = min(epoch)-5 & xmax = max(epoch)+5
;xmin = -75 & xmax = 75
plot, fix(epoch), omc,psym=3, yrange=[ymin,ymax],xrange=[xmin,xmax],$
  xtitle='Epoch',ytitle='O-C (min)';, /xstyle;,xticks=4
oploterr, epoch, omc, err, 3
oplot, [-9d9,9d9],[0,0],linestyle=1

syms = [0,3,8,5,0,3,8,5]
fill = [1,1,1,1,0,0,0,0]
nsyms = n_elements(syms)

;; overplot each telescope in different colors
if n_elements(telescope) ge 1 then begin
   sorted = sort(telescope)
   tnames = telescope[sorted[uniq(telescope[sorted])]]
   for i=0, n_elements(tnames)-1 do begin
      observed = where(telescope eq tnames[i])
      if observed[0] ne -1 then begin
         plotsym, syms[i mod nsyms], color=colors[i mod ncolors],fill=fill[i mod nsyms]
         oplot, epoch[observed],omc[observed],psym=8
         xsize = (!x.crange[1] - !x.crange[0])
         ysize = (!y.crange[1] - !y.crange[0])
         xyouts, !x.crange[0] + xlegend*xsize,!y.crange[0]+(ylegend - i*charsizelegend)*ysize, $
                 tnames[i],color=colors[i mod ncolors],charsize=charsize
         oplot, [!x.crange[0]+xlegend*xsize-xsize/20],$
                [!y.crange[0]+(ylegend - (i-0.25)*charsizelegend)*ysize],psym=8
      endif
   endfor
endif

;; print the best-fit ephemeris, uncertainties, and chi^2
print, 'T0 = '+string(coeffs[0],format='(f14.6)')+' +/- '+strtrim(sigma[0],2)
print, 'P  = '+strtrim(coeffs[1],2) + ' +/- ' + strtrim(sigma[1],2)
print, 'chi^2 = '+strtrim(total(((omc)/err)^2),2)
print, 'dof = '+strtrim((n_elements(time)-2),2)
print, 'chi^2/dof = ' + strtrim(total(((omc)/err)^2)/(n_elements(time)-2),2)
print, ''

;; print the epochs, times, errors, etc.
;print, 'Epoch      T_C      ERR    O-C (O-C)/err'
;exofast_forprint, epoch, time, err*60, omc*60d0, omc/err,/t, $
;  format='(i4,x,f14.6,x,i4,x,f7.2,x,f5.2)'

;; close the postscript device
if keyword_set(ps) or n_elements(epsname) ne 0 then begin
   !p.font=-1
   !p.multi=0
   device, /close
   device, encapsulated=0
   set_plot, mydevice
endif

end
