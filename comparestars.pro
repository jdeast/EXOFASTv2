;+
; NAME:
;   COMPARESTARS
;
; PURPOSE: 
;   Compares the stellar parameters of several different fits using
;   the IDL save files as output by exofastv2.pro.
;
; CALLING SEQUENCE:
;   comparestars, 'savpath' [,TAGS=,PSNAME=]
;
; INPUTS:
;   SAVPATH - A string resolved by file_search to specify the names of
;             all the idl save files to compare.
;   TAGS    - A string array specifying the tagnames (parameters) to
;             compare. Default is 
;             ['mstar','rstar','rhostar','logg','teff','feh','initfeh','age','lstar','eep']
;   PSNAME  - A string specifying the name of the postscript file to
;             output. If not specified, the plots are made to the screen.
;-
pro comparestars, savpath, tags=tags, psname=psname

;; these are the (stellar) tags that will be compared between each file
if n_elements(tags) eq 0 then $
   tags = ['mstar','rstar','rhostar','logg','teff','feh','initfeh','age','lstar','eep']

files = file_search(savpath,count=nfiles)

if n_elements(psname) ne 0 then begin
   mydevice = !d.name
   set_plot, 'PS'
   loadct, 39, /silent
   !p.font=0
   aspect_ratio=1.5
   xsize=9
   ysize=xsize/aspect_ratio
   device,  filename=psname,/color,bits=24
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
   colors = [black,red,green,blue,orange,cyan,purple,redorange,yellowgreen,$
             darkpurple,lightgreen,greenblue,darkblue,bluegreen,lightblue,yellow]
   legendsize = 1
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

;; create an array of structures from the idl save files
structures = ptrarr(nfiles,/allocate_heap)
labels = strarr(nfiles)
for j=0L, nfiles-1 do begin
   restore, files[j]
   *structures[j] = mcmcss
   labels[j] = (strsplit(files[j],'.',/extract))[1]
endfor

;; make a histogram of the parameter for each save file, one parameter per page
;; for each stellar parameter....
for i=0, n_elements(tags)-1 do begin

   ;; find the indicies in the structure that corresponds to the parameter
   tagmatch = lonarr(nfiles)
   xmax = -!values.d_infinity
   xmin = !values.d_infinity
   for j=0, nfiles-1 do begin
      tagmatch[j] = (where(tag_names((*structures[j]).star) eq strupcase(tags[i])))[0]
      if tagmatch[j] eq -1 then begin
         print, 'No match for tag ' + tags[i]
         continue
      endif
      ;; and find the range of the parameter
      burnndx = (*structures[j]).burnndx
      nchains = (*structures[j]).nchains
      nsteps = (*structures[j]).nsteps/nchains
    
      chi2 = reform((*(*structures[j]).chi2),nsteps,nchains)
      burnndx = getburnndx(chi2,goodchains=goodchains)
      par = (reform((*structures[j]).star.(tagmatch[j]).value,nsteps,nchains))[burnndx:*,goodchains]

      thismax = max(par,min=thismin)
      if thismax gt xmax then xmax = thismax
      if thismin lt xmin then xmin = thismin
   endfor

   ;; make a histogram of the parameter for each save file
   hists = dblarr(nfiles,100)
   for j=0, nfiles-1 do begin
      burnndx = (*structures[j]).burnndx

      nchains = (*structures[j]).nchains
      nsteps = (*structures[j]).nsteps/nchains
      chi2 = reform((*(*structures[j]).chi2),nsteps,nchains)
      burnndx = getburnndx(chi2,goodchains=goodchains)
      
      par = (reform((*structures[j]).star.(tagmatch[j]).value,nsteps,nchains))[burnndx:*,goodchains]

      hist = histogram(par,nbin=100,locations=x,max=xmax,min=xmin)
      hists[j,*] = hist/total(hist)
      ;; zero out singular arrays (e.g., age when age isn't fit)
      bad = where(hists[j,*] eq 1d0)
      if bad[0] ne -1 then hists[j,*] = 0d0
   endfor

   ;; plot all the histograms (with a legend)
   ymax=max(hists)
   plot, [0], [0], xrange=[xmin,xmax], yrange=[0,ymax], xtitle=exofast_textoidl((*structures[0]).star.(tagmatch[0]).latex), ytitle='Normalized Probability'
   for j=0L, nfiles-1 do oplot, x, hists[j,*], psym=10, color=colors[j mod ncolors], thick=thick
   exofast_legend, labels,/top,/right,color=colors[lindgen(nfiles) mod ncolors],linestyle=lonarr(nfiles),charsize=legendsize

endfor

;; make a covariance plot of each parameter against all others for each save file, one plot
;; for each stellar parameter (x axis)....
for i=0, n_elements(tags)-1 do begin

   ;; find the indicies in the structure that corresponds to the parameter
   tagmatchx = lonarr(nfiles)
   xmax = -!values.d_infinity
   xmin = !values.d_infinity
   for j=0, nfiles-1 do begin
      tagmatchx[j] = (where(tag_names((*structures[j]).star) eq strupcase(tags[i])))[0]
      if tagmatchx[j] eq -1 then begin
         print, 'No match for tag ' + tags[i]
         continue
      endif
      ;; and find the range of the paramete
      burnndx = (*structures[j]).burnndx

      nchains = (*structures[j]).nchains
      nsteps = (*structures[j]).nsteps/nchains
      chi2 = reform((*(*structures[j]).chi2),nsteps,nchains)
      burnndx = getburnndx(chi2,goodchains=goodchains)
      
      par = (reform((*structures[j]).star.(tagmatchx[j]).value,nsteps,nchains))[burnndx:*,goodchains]

      thismax = max(par,min=thismin)
      if thismax gt xmax then xmax = thismax
      if thismin lt xmin then xmin = thismin
   endfor
   
   ;; for each other stellar parameter (y axis)....
   for k=i+1, n_elements(tags)-1 do begin

      ;; find the indicies in the structure that corresponds to the parameter
      tagmatchy = lonarr(nfiles)
      ymax = -!values.d_infinity
      ymin = !values.d_infinity
      for j=0, nfiles-1 do begin
         tagmatchy[j] = (where(tag_names((*structures[j]).star) eq strupcase(tags[k])))[0]
         if tagmatchy[j] eq -1 then begin
            print, 'No match for tag ' + tags[k]
            continue
         endif
         ;; find the range of this variable
         burnndx = (*structures[j]).burnndx
 
         nchains = (*structures[j]).nchains
         nsteps = (*structures[j]).nsteps/nchains
         chi2 = reform((*(*structures[j]).chi2),nsteps,nchains)
         burnndx = getburnndx(chi2,goodchains=goodchains)

         par = (reform((*structures[j]).star.(tagmatchy[j]).value,nsteps,nchains))[burnndx:*,goodchains]
      
         thismax = max(par,min=thismin)
         if thismax gt ymax then ymax = thismax
         if thismin lt ymin then ymin = thismin
      endfor

      ;; draw the shell of the plot with the right ranges and axis labels
      plot, [0],[0],xrange=[xmin,xmax],yrange=[ymin,ymax], $
            xtitle=exofast_textoidl((*structures[0]).star.(tagmatchx[0]).latex),$
            ytitle=exofast_textoidl((*structures[0]).star.(tagmatchy[0]).latex),$
            /xstyle,/ystyle

      ;; overplot get the covariance contours for each save file
      for j=0, nfiles-1 do begin
         burnndx = (*structures[j]).burnndx
         nchains = (*structures[j]).nchains
         nsteps = (*structures[j]).nsteps/nchains

         chi2 = reform((*(*structures[j]).chi2),nsteps,nchains)
         burnndx = getburnndx(chi2,goodchains=goodchains)

         x = (reform((*structures[j]).star.(tagmatchx[j]).value,nsteps,nchains))[burnndx:*,goodchains]
         y = (reform((*structures[j]).star.(tagmatchy[j]).value,nsteps,nchains))[burnndx:*,goodchains]

         if min(x) ne max(x) and min(y) ne max(y) then begin
            exofast_errell,x,y,xpath=xpath,ypath=ypath,$
                           prob=probs,nxbin=100, nybin=100,$
                           xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,logname=logname
            if n_elements(xpath) gt 0 and n_elements(ypath) eq n_elements(xpath) then $
               oplot, xpath, ypath, color=colors[j mod ncolors], thick=thick
         endif
      endfor

      ;; now plot the legend
      exofast_legend, labels,/top,/right,color=colors[lindgen(nfiles) mod ncolors],$
                      linestyle=lonarr(nfiles),charsize=legendsize


   endfor ;; each parameter on the Y axis
endfor ;; each parameter on the X axis

if keyword_set(ps) then begin
   !p.font=-1
   !p.multi=0
   device, /close
   set_plot, mydevice
;   spawn, 'gv ' + psname + ' &'
endif

end
