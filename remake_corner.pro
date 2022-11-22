;+
; NAME:
;   REMAKE_CORNER
;
; PURPOSE: 
;   Creates a corner plot of several different fits overlaid on top
;   of each other using the IDL save files as output by exofastv2.pro.
;
; CALLING SEQUENCE:
;   remake_corner, 'savpath' [,TAGS=,LATEXNAMES=,LEGENDTXT=,PSNAME=]
;
; INPUTS:
;  SAVPATH    - A string resolved by file_search or an array of
;               filenames specifying the names of all the idl save
;               files to compare.  
;               See NSAMPLE if files are large
;  TAGS       - A string array specifying the tagnames (parameters) to
;               compare. Default is
;               ['mstar','rstar','rhostar','logg','teff','feh','initfeh','age','lstar','eep']
;  LATEXNAMES - A string array specifying the textoidl axis labels
;               that corresponds to each of the tags above. If not the
;               same length as TAGS, it will use the structure's latex
;               property, but that doesn't always look pretty in
;               textoidl. Empty strings will use the the structure's
;               latex property.
;  PSNAME     - A string specifying the name of the postscript file to
;               output. If not specified, the plots are made to the
;               screen.
;  LEGENDTXT  - A string array corresponding to each save file that
;               specifies the legend txt. Default is the filenames.
;  NSAMPLE    - A scalar that specifies how many (random) samples to
;               take from each parameter. Default uses all. Useful for
;               memory management if plotting data from many large files
;  POSITION   - A 2-element array specifying the normalized X and Y
;               coordinates of the top left corner of the legend
;  POSTERIORS - When this keyword is set, text files of the posteriors
;               will be output with the same name as the input files,
;               replacing 'idl' with 'txt'.
;
;-
pro remake_corner, savpath, tags=tags, latexnames=latexnames, psname=psname, $
                   legendtxt=legendtxt, nsample=nsample, position=position, $
                   posteriors=posteriors, nxbin=nxbin, nybin=nybin, optimizeorder=optimizeorder

if n_elements(position) ne 2 then position = [0.533d0,0.903d0]

;; these are the (stellar) tags that will be compared between each file
if n_elements(tags) eq 0 then $
   tags = ['mstar','rstar','rhostar','logg','teff','feh','initfeh','age','lstar','eep']

npars = n_elements(tags)

if n_elements(latexnames) ne npars then latexnames = strarr(npars)

if n_elements(savpath) gt 1 then begin
   files = savpath
   nfiles = n_elements(savpath)
endif else files = file_search(savpath,count=nfiles)

if n_elements(psname) ne 0 then begin
   mydevice = !d.name
   set_plot, 'PS'
   loadct, 39, /silent
   !p.font=0
   aspect_ratio=1.5
   xsize=9
   ysize=xsize/aspect_ratio
   device,  filename=psname,/color,bits=24
   if npars eq 1 then device, xsize=xsize, ysize=ysize
   red = 254
   green = 159
   purple = 32
   cyan = 104
   orange = 201
   yellow = 207
   bluegreen = 111
   lightgreen = 191

   black = 0         ; 0
   darkpurple = 16   ; 9
   purple = 32       ; 6
   darkblue = 48     ; 12 
   blue = 64         ; 3
   lightblue = 80    ; 14
   cyan = 96         ; 5
   bluegreen = 112   ; 13
   greenblue = 128   ; 11
   green = 144       ; 2
   lightgreen = 160  ; 10
   yellowgreen = 176 ; 8 
   yellow = 192      ; 15
   orange = 208      ; 4
   redorange = 224   ; 7
   red = 254         ; 1
   colors = [black,red,blue, green,orange,cyan,purple,redorange,yellowgreen,$
             darkpurple,lightgreen,greenblue,darkblue,bluegreen,lightblue,yellow]
   legendsize = 0.5
   thick=2
endif else begin
   red = '0000ff'x
   green = '00ff00'x
   blue = 'ff0000'x
   purple = 'ff00ff'x
   cyan = 'ffff00'x
   white = 'ffffff'x
   legendsize = 1.5
   colors = [red,green,blue,purple,cyan,white]
endelse
ncolors = n_elements(colors)

;; extract all info we'll need from each save file
;; this can eat up a lot of memory!
labels = strarr(nfiles)
pars = ptrarr(nfiles,npars,/allocate_heap)
extrema = dblarr(2,nfiles,npars) + !values.d_nan

for j=0L, nfiles-1 do begin
   restore, files[j]

   labels[j] = (strsplit(files[j],'.',/extract))[1]
   
   burnndx = mcmcss.burnndx
   nchains = mcmcss.nchains
   nsteps = mcmcss.nsteps/nchains
   goodchains = (*mcmcss.goodchains)
   ngoodsteps = n_elements(goodchains)*(nsteps-burnndx)

   ;; save memory by randomly sampling
   if n_elements(nsample) ne 0 then sample = long(randomu(seed,nsample)*ngoodsteps) $
   else sample = lindgen(ngoodsteps) ;; or not

   if keyword_set(posteriors) then pars2txt = dblarr(npars,ngoodsteps)

   ;; extract the relevant values from the structure
   for i=0L, npars-1 do begin

      parstr = find_by_tag(mcmcss, tags[i])
      if (size(parstr))[2] ne 8 then continue

      par = ((reform(parstr.value,nsteps,nchains))[burnndx:*,goodchains])[sample]
      xmin = min(par,max=xmax)
      extrema[0,j,i] = xmin
      extrema[1,j,i] = xmax
      if latexnames[i] eq '' then latexnames[i] = parstr.latex

      *pars[j,i] = par
      if keyword_set(posteriors) then pars2txt[i,*] = par

   endfor
   undefine, mcmcss ;; free up memory before we load the next one

   if keyword_set(posteriors) then begin
      openw, lun, file_basename(files[j],'.idl') + '.txt',/get_lun
      printf, lun, tags, format='("#"' + strtrim( npars,2) + '(a27,x))'
      printf, lun, pars2txt, format='(' + strtrim( npars,2) + '(f27.20,x))'
      free_lun, lun
      undefine, pars2txt
   endif
endfor

if keyword_set(optimizeorder) then begin
   area = dblarr(nfiles,npars,npars)
   area2 = dblarr(nfiles,npars,npars)
   if n_elements(nxbin) eq 0 then nxbin = 100
   if n_elements(nybin) eq 0 then nybin = nxbin
   
   for j=0L, nfiles-1 do begin   
      for i=0L, npars-1 do begin
         xpar = *pars[j,i]
         xmin = min(extrema[0,*,i])
         xmax = max(extrema[1,*,i])
         for k=0L, npars-1 do begin
            ypar = *pars[j,k]
            ymin = min(extrema[0,*,k])
            ymax = max(extrema[1,*,k])
            xbinsz = double(xmax-xmin)/(nxbin-1)
            ybinsz = double(ymax-ymin)/(nybin-1)
            hist2d = hist_2d(xpar,ypar,$
                             bin1=xbinsz,min1=xmin,max1=xmax,$
                             bin2=ybinsz,min2=ymin,max2=ymax)
            hist2d = hist2d/total(hist2d)
            sorted = sort(hist2d)
            thresh = hist2d[sorted[n_elements(hist2d)*0.05]]
;            area[j,i,k] = total(1d0/hist2d^2,/nan)
            area[j,i,k] = n_elements(where(hist2d gt thresh))
         endfor
      endfor
   endfor
endif
;stop

if n_elements(legendtxt) ne nfiles then legendtxt = labels

;; prepare the corner plot
!p.multi=0
npars = n_elements(tags)
if npars gt 1 then begin
   charsize = 4d0/npars
   exofast_multiplot, [npars,npars], /square, mXtitle='', mYtitle='', /rowmajor
endif

;; for each parameter (x axis)....
for i=0L, n_elements(tags)-1 do begin
   xmin = min(extrema[0,*,i],/nan)
   xmax = max(extrema[1,*,i],/nan)
   
   ;; for each other parameter (y axis)....
   for k=0L, n_elements(tags)-1 do begin

      if i eq 0 then ytitle='!3' + exofast_textoidl(latexnames[k]) else ytitle = ''
      if k eq npars-1 then xtitle='!3' + exofast_textoidl(latexnames[i]) else xtitle = ''
      if npars eq 1 then begin
         ytitle='Normalized Prob'
         charsize=1
      endif else begin
         xminor=1
         yminor=1
         xstyle=1
         xtickname=replicate(" ",60)
         ytickname=replicate(" ",60)
      endelse

      ;; if they're the same parameter, make a histogram
      if k lt i+1 then begin
         if k eq i then begin

            ;; make a histogram of the parameter for each save file
            hists = dblarr(nfiles,100)
            for j=0, nfiles-1 do begin

               par = *pars[j,k]
               if n_elements(par) eq 0 then continue

;xmin = 1.6
;xmax = 1.9
               hist = histogram(par,nbin=100,locations=x,min=xmin, max=xmax)
               hists[j,*] = hist/total(hist)
               ;; zero out singular arrays (e.g., age when age isn't fit)
               bad = where(hists[j,*] eq 1d0)
               if bad[0] ne -1 then hists[j,*] = 0d0
            endfor
            
            ;; plot all the histograms
            ymax=max(hists)
            plot, [0], [0], xrange=[xmin,xmax], yrange=[0,ymax], xtitle=xtitle, ytitle=ytitle,$
                  xminor=xminor, yminor=yminor,xtickname=xtickname,ytickname=ytickname,$
                  charsize=charsize,xstyle=xstyle
            for j=0L, nfiles-1 do oplot, x, hists[j,*], psym=10, color=colors[j mod ncolors], thick=thick

         endif
         exofast_multiplot ;; next plot
         continue ;; skip other half of corner plot
      endif

      ;; make the contour of the covariance
      ymin = min(extrema[0,*,k],/nan)
      ymax = max(extrema[1,*,k],/nan)

      ;; draw the shell of the plot with the right ranges and axis labels
      plot, [0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax], $
            thick=thick,charsize=charsize, $
            /ystyle,/xstyle,xminor=1, yminor=1,$ 
            xtickname=replicate(" ",60),ytickname=replicate(" ",60), xtitle=xtitle, ytitle=ytitle
      
      ;; plot them from largest are to smallest, so they don't
      ;; cover each other up
      if keyword_set(optimizeorder) then begin
         order = reverse(sort(area[*,i,k]))
      endif else order = lindgen(nfiles)

      ;; overplot get the covariance contours for each save file
      for j=0, nfiles-1 do begin

         x = *pars[order[j],i]
         y = *pars[order[j],k]

         if min(x) ne max(x) and min(y) ne max(y) then begin
            exofast_errell,x,y,xpath=xpath,ypath=ypath,$
                           prob=probs,nxbin=nxbin, nybin=nybin,$
                           xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,logname=logname
            if n_elements(xpath) gt 0 and n_elements(ypath) eq n_elements(xpath) then $
               oplot, xpath, ypath, color=colors[order[j] mod ncolors], thick=thick
         endif
      endfor
      exofast_multiplot ;; next plot

   endfor ;; each parameter on the Y axis
endfor ;; each parameter on the X axis

;; draw the legend
exofast_legend, legendtxt,color=colors[lindgen(nfiles) mod ncolors],$
                linestyle=lonarr(nfiles),charsize=legendsize,pos=position,/norm

;; reset back to default plotting parameters
exofast_multiplot, /reset
exofast_multiplot, /default 

if n_elements(psname) ne 0 then begin
   !p.font=-1
   !p.multi=0
   device, /close
   set_plot, mydevice
;   spawn, 'evince ' + psname + ' &'
;   spawn, 'convert -flatten -density 300x300 ' + psname + ' ' + file_basename(psname,'.ps') + '.png'
endif

end
