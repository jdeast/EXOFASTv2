pro mksummaryframe, ss0, ndx, base=base, idlfile=idlfile, printthisval=printthisval,redo=redo

  COMMON chi2_block, ss
  
  if n_elements(idlfile) ne 0 then begin
     ;; this is the (1 element) structure used by exofast_chi2v2
     restore, idlfile
     ss = mcmcss2ss(mcmcss)
;     junk = step2pars(ss)
;     stop
;     junk = pars2step(ss)
     if n_elements(ndx) eq 0 then minchi2 = min(*mcmcss.chi2,ndx)
     pars = str2pars(mcmcss,ndx=ndx)
  endif else begin
     ss = ss0
     pars = str2pars(ss,ndx=ndx)
  endelse
  
  chi2 = exofast_chi2v2(pars,psname=base)
  
  evolutionarymodel = 'mist'
  summaryfile = base + '.exofastv2.png'

  extractpgno, base+'.transit.ps', stacked_primary=stackedpage, primary=transitpage
  
  ;; convert multi-page ps files to pdfs, crop, convert to png
  if file_test(base+'.transit.ps') then begin
     spawn, 'ps2pdf ' + base + '.transit.ps '  + base + '.transit.pdf'
     spawn, 'grep "PageBoundingBox:" ' + base + '.transit.ps', boundingboxes
  endif
  
  if file_test(base+'.rv.ps') then $
     spawn, 'ps2pdf ' + base + '.rv.ps '  + base + '.rv.pdf'

  for i=0L, ss.nplanets-1 do begin
     
     ;; *********** phased transit plot *********
     pageno = strtrim(transitpage[i],2)

     ;; can these be combined into one step?
     ;; extract the relevant page into single page pdf file
     spawn, 'gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -DSAFER -dFirstPage=' +$
            pageno + ' -dLastPage=' + pageno + $
            ' -sOutputFile='+base+'.transit.'+strtrim(i,2)+'.pdf ' +$
            base + '.transit.pdf', output

     ;; crop pdf file
     spawn, 'convert -density 600x600 ' +$
            base + '.transit.' + strtrim(i,2) + '.pdf '+$
            '-crop 2600x1650+350+2325 '+$
            '+repage ' + base + '.transit.' + strtrim(i,2) + '.cropped.pdf', output

     ;; convert to png file
     spawn, 'convert -flatten -density 300x300 ' + $
            base + '.transit.' + strtrim(i,2) + '.cropped.pdf ' +$
            base + '.transit.' + strtrim(i,2)+ '.png', output
     ;; *****************************************

     ;; ********** stacked transit plot *********
     pageno = strtrim(stackedpage[i],2)
     oversamp = 8
     arr = long((strsplit(boundingboxes[pageno-1],/extract))[1:4])
     width = strtrim((arr[2]-arr[0])*oversamp,2)
     height = strtrim((arr[3]-arr[1])*oversamp,2)
     x1 = strtrim(arr[0]*oversamp,2)
     y1 = strtrim((836-arr[3])*oversamp,2)
     density = strtrim(72L*oversamp,2)

     ;; extract the relevant page into a single page pdf file
     spawn, 'gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -DSAFER -dFirstPage=' +$
            pageno + ' -dLastPage=' + pageno +$
            ' -sOutputFile='+base+'.transit.stacked.'+strtrim(i,2)+'.pdf '+$
            base + '.transit.pdf', output

     ;; crop the pdf file
     spawn, 'convert -density ' + density + 'x' + density + ' ' +$
            base + '.transit.stacked.' + strtrim(i,2) + '.pdf '+$
            '-crop ' + width + 'x' + height + '+' + x1 + '+' + y1 +$
            ' +repage '+base+'.transit.stacked.'+strtrim(i,2)+'.cropped.pdf', output

     ;; convert to png
     spawn, 'convert -flatten -density 300x300 ' +$
            base + '.transit.stacked.' + strtrim(i,2) + '.cropped.pdf ' +$
            base + '.transit.stacked.' + strtrim(i,2)+ '.png', output
     ;; *****************************************

     ;; *********** phased rv plot *********
     if file_test(base+'.rv.ps') then begin

        pageno = strtrim(i+1,2)

        ;; can these be combined into one step?
        ;; extract the relevant page into single page pdf file
        spawn, 'gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -DSAFER -dFirstPage=' +$
               pageno + ' -dLastPage=' + pageno + $
               ' -sOutputFile='+base+'.rv.'+strtrim(i,2)+'.pdf ' +$
               base + '.rv.pdf', output
        
        ;; crop pdf file
        spawn, 'convert -density 600x600 ' +$
               base + '.rv.' + strtrim(i,2) + '.pdf '+$
               '-crop 2300x1650+600+2325 '+$
               '+repage ' + base + '.rv.' + strtrim(i,2) + '.cropped.pdf', output
        
        ;; convert to png file
        spawn, 'convert -flatten -density 300x300 ' + $
               base + '.rv.' + strtrim(i,2) + '.cropped.pdf ' +$
               base + '.rv.' + strtrim(i,2)+ '.png', output
        ;; *****************************************
     endif
        
  endfor
  
  ;; extract best parameters to print
  annotation = strarr(ss.nplanets+ss.nstars)

  ;; stellar parameters to include
  varnames = ['mstar','rstar','teff','feh','logg','errscale','Av','age','parallax','dilute']
  for i=0L, ss.nstars-1 do begin
     for j=0L, n_elements(varnames)-1 do begin
        parname = varnames[j] + '_' + strtrim(i,2)

        map = findvariable(ss,parname,count=0)

        ;; extract the parameter using the map
        if map[4] ne -1 then begin
           parameter = (*ss.(map[0])[map[1]].(map[2])[map[3]]).(map[4])[map[5]]
           mcmcpar = (*mcmcss.(map[0])[map[1]].(map[2])[map[3]]).(map[4])[map[5]]
        endif else if map[2] ne -1 then begin 
           parameter = ss.(map[0])[map[1]].(map[2])[map[3]]
           mcmcpar = mcmcss.(map[0])[map[1]].(map[2])[map[3]]
        endif else if map[0] ne -1 then begin
           parameter = ss.(map[0])[map[1]]
           mcmcpar = mcmcss.(map[0])[map[1]]
        endif

        ;; median + 68% CI
;        value = [parameter.medvalue,parameter.upper,parameter.lower]
;        annotation[i] += string(parname,value, format='(a,a,"+",a,"-",a,"\n")')

        ;; value at this step
        value = parameter.value;[ndx]

        if keyword_set(printthisval) then begin
           annotation[i] += string(parname, value, format='(a,a,"\n")')
        endif else begin
;           annotation[i] += string(parname, value, format='(a,a)')
           annotation[i] += parname + ' ' + mcmcpar.medvalue + ' +' + mcmcpar.upper + ' -' + mcmcpar.lower + '\n'
        endelse
           
     endfor
  endfor

  ;; planetary parameters to include
  varnames = ['rp','Period','b','e','omegadeg','ideg']
  for i=0L, ss.nplanets-1 do begin
     for j=0L, n_elements(varnames)-1 do begin
        parname = varnames[j] + '_' + strtrim(i,2)
        
        map = findvariable(ss,parname,count=0)

        ;; extract the parameter using the map
        if map[4] ne -1 then begin
           parameter = (*ss.(map[0])[map[1]].(map[2])[map[3]]).(map[4])[map[5]]
           mcmcpar =  (*mcmcss.(map[0])[map[1]].(map[2])[map[3]]).(map[4])[map[5]]
        endif else if map[2] ne -1 then begin 
           parameter = ss.(map[0])[map[1]].(map[2])[map[3]]
           mcmcpar = mcmcss.(map[0])[map[1]].(map[2])[map[3]]
        endif else if map[0] ne -1 then begin
           parameter = ss.(map[0])[map[1]]
           mcmcpar = mcmcss.(map[0])[map[1]]
        endif
        
        ;; median + 68% CI
;        value = [parameter.medvalue,parameter.upper,parameter.lower]
;        annotation[ss.nstars+i] += string(parname,value, format='(a,a,"+",a,"-",a,"\n")')

        ;; value at this step
        value = parameter.value;[ndx]
        if keyword_set(printthisval) then begin
           annotation[ss.nstars+i] += string(parname, value, format='(a,a,"\n")')
        endif else begin
           ;annotation[ss.nstars+i] += string(parname, value, format='(a,a)')
           annotation[ss.nstars+i] += parname + ' ' + mcmcpar.medvalue + ' +' + mcmcpar.upper + ' -' + mcmcpar.lower + '\n'
        endelse
     endfor
  endfor

  ;; create a blank white canvas with room for all planets
  nx = 3                               ;; number of columns
  ny = 1+ss.nplanets+ss.nstars ;; number of rows
  xsize = 1000                         ;; xsize of the canvas
  ysize = xsize*8.5d0/11d0*0.5/2d0*ny  ;; ysize of the canvas
  spawn, 'convert -size ' + strtrim(xsize*0.95,2) + 'x' + strtrim(ysize,2)+ $
         ' xc:white ' + summaryfile, output

  ;; make a row for each planet showing 
  ;;   column 1: the parameters of the model
  ;;   column 2: phased transit model
  ;;   column 3: phased RV model
  for i=0L, ss.nplanets-1 do begin

     ;; add the important planet parameters as text to the left side 
     ;; in line with each planet
     spawn, 'convert -gravity SouthWest ' + summaryfile + $
            ' -font NimbusSans-Regular' +$
            ' -annotate +100+' + strtrim(75+(i+2)*190,2) + " '"  + $
            annotation[i+1] + "' " + summaryfile, output 

     ;; add the phased transit to the canvas
     phasedpngname = base + '.transit.'+strtrim(i,2) + '.png'   
     spawn, 'composite -geometry ' + $
            strtrim(round(xsize/nx)*(1),2)+'x'+$ ;; x width (pixels)
            strtrim(round(ysize/ny)*(1),2)+'+'+$ ;; y width (pixels)
            strtrim(round(xsize/nx)*(1),2)+'+'+$ ;; x offset (pixels)
            strtrim(round(ysize/ny)*(i),2)+' '+$ ;; y offset (stack planets in Y)
            phasedpngname + ' ' + $              ;; file to overlay
            summaryfile + ' '  + summaryfile,$   ;; create composite image
            output
     
     ;; add the phased RV plot to the canvas
     rvpngname = base + '.rv.'+strtrim(i,2) + '.png'
     if file_test(rvpngname) then begin
        cmd = 'composite -geometry ' + $
              strtrim(round(xsize/nx)*(1),2)+'x'+$ ;; x width (pixels)
              strtrim(round(ysize/ny)*(1),2)+'+'+$ ;; y width (pixels)
              strtrim(round(xsize/nx)*(2),2)+'+'+$ ;; x offset (pixels)
            strtrim(round(ysize/ny)*(i),2)+' '+$ ;; y offset (stack planets in Y)
              rvpngname + ' ' + $                ;; file to overlay
              summaryfile + ' '  + summaryfile   ;; create composite image
        spawn, cmd, output
     endif else begin
        ;; add the stacked (per LC file) transits to the canvas
        stackedpngname = base + '.transit.stacked.'+strtrim(i,2) + '.png'
        spawn, 'composite -geometry ' + $
               strtrim(round(xsize/nx)*(1),2)+'x'+$ ;; x width (pixels)
               strtrim(round(ysize/ny)*(1),2)+'+'+$ ;; y width (pixels)
               strtrim(round(xsize/nx)*(2),2)+'+'+$ ;; x offset (pixels)
               strtrim(round(ysize/ny)*(i),2)+' '+$ ;; y offset (stack planets in Y)
               stackedpngname + ' ' + $             ;; file to overlay
               summaryfile + ' '  + summaryfile     ;; create composite image
     endelse
     
  endfor

  for i=0L, ss.nstars-1 do begin
     ;; overlay stellar parameter text
     cmd = 'convert -gravity SouthWest ' + summaryfile +$
            ' -font NimbusSans-Regular' +$
            ' -annotate +100+230 ' + "'" + annotation[i] + "' " +$
            summaryfile
     spawn, cmd, output

     ;; overlay the sed plot 
     spawn, 'convert -flatten -density 300x300 ' + base + '.sed.eps ' + base + '.sed.png' 
     cmd = 'composite -geometry ' + $
            strtrim(round(xsize/nx)*1,2) + 'x' +$
            strtrim(round(ysize/ny)*1,2) + '+' +$
            strtrim(round(xsize/nx)*1,2) + '+' +$
            strtrim(round(ysize/ny*(ss.nplanets+i)),2)+' '+$
            base + '.sed.png ' +$
            summaryfile + ' '  + summaryfile
;     print, cmd
     spawn, cmd, output

     ;; overlay the evolutionary model
     emfilebase = base + '.' + evolutionarymodel + '.' + string(i,format='(i03)')
     spawn, 'convert -flatten -density 300x300 ' + emfilebase + '.eps ' + emfilebase + '.png' 
     cmd = 'composite -geometry ' + $
           strtrim(round(xsize/nx*0.68),2) + 'x' +$
           strtrim(round(ysize/ny*1.00),2) + '+' +$
           strtrim(round(xsize/nx*2.00),2) + '+' +$
           strtrim(round(ysize/ny*(ss.nplanets+i)),2) +' ' +$
           emfilebase + '.png ' +$
           summaryfile + ' '  + summaryfile
;     print, cmd
     spawn, cmd, output

  endfor

  
  if not file_test(base + '.planet_corner.png') or keyword_set(redo) then begin
     remake_corner, idlfile, tags=['mstar','rstar','teff','av','age'     ],$
                    sample=sample, psname=base+'.star_corner.ps',extrema=starextrema
     spawn, 'convert -flatten -density 300x300 ' + base + '.star_corner.ps ' + base + '.star_corner.png'
  endif
  cmd = 'composite -geometry ' +$
        strtrim(round(xsize/nx*1.00),2) + 'x' +$
        strtrim(round(ysize/ny*1.00),2) + '+' +$
        strtrim(round(xsize/nx*0.25),2) + '+' +$
        strtrim(round(ysize/ny*(ss.nplanets+ss.nstars)),2)+' '+$
        base + '.star_corner.png '+$ 
        summaryfile + ' '  + summaryfile
  spawn, cmd, output

  ;; make a corner plot for the planet parameters
  if ss.nplanets ge 1 then begin
     if not file_test(base + '.planet_corner.png') or keyword_set(redo) then begin
        remake_corner, idlfile, tags=['mp'   ,   'rp', 'b'  ,'e' ,'omegadeg'],$
                       sample=sample, psname=base+'.planet_corner.ps',extrema=planetextrema
        spawn, 'convert -flatten -density 300x300 ' + base + '.planet_corner.ps ' + base + '.planet_corner.png'
     endif
     cmd = 'composite -geometry ' +$
           strtrim(round(xsize/nx*1.00),2) + 'x' +$
           strtrim(round(ysize/ny*1.00),2) + '+' +$
           strtrim(round(xsize/nx*2.00),2) + '+' +$
           strtrim(round(ysize/ny*(ss.nplanets+ss.nstars)),2)+' '+$
           base + '.planet_corner.png '+$ 
           summaryfile + ' '  + summaryfile
     spawn, cmd, output
  endif
     
  ;; Likelihood chain plot
  chain_png_name = base + '.chi2chain.png'
  burnndx = mcmcss.burnndx
  goodchains = *mcmcss.goodchains
  nchains = mcmcss.nchains
  rgb = colortable(39,ncolors=nchains+1)
  white = 256L^2*255L + 256L*255 + 255
  black = 0L

  xsizechain=500
  ysizechain=200
  set_plot, 'Z',/copy
  device, set_resolution=[xsizechain,ysizechain], set_pixel_depth=24, set_font='Times',/tt_font
  chi2 = reform(*mcmcss.chi2,mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
  loglike = -chi2/2d0

  ymin = min(loglike,max=ymax)
  if ymin lt 0 then sign = ' + ' else sign = ' - '
  latex = textoidl('\Delta log(Like)') ; + sign + strtrim(abs(ymax),2)                                                                            
  yminplot = min(loglike[burnndx:-1,goodchains],max=ymax)
  yminplot = -40+ymax
  ymaxplot = 0d0+ymax
  nsteps = mcmcss.nsteps/mcmcss.nchains
  plot, [0],[0], xrange=[0,nsteps-1],yrange=[yminplot,ymaxplot]-ymax, $
        ytitle=latex,xtitle='Chain Link number',/xs,background=white,color=black,/ys
  for j=0, n_elements(goodchains)-1 do begin
     color = rgb[goodchains[j],0]*256L^2 + rgb[goodchains[j],1]*256L + rgb[goodchains[j],2]
     oplot, loglike[0:nsteps-1,goodchains[j]]-ymax,color=color
  endfor
  oplot, [burnndx,burnndx],[-9d9,9d9],color=black
  write_png, chain_png_name, tvrd(/true)
  cmd = 'composite -geometry ' +$
        strtrim(round(xsize/nx*1.00),2) + 'x' +$
        strtrim(round(ysize/ny*1.00),2) + '+' +$
        strtrim(round(xsize/nx*1.00),2) + '+' +$
        strtrim(round(ysize/ny*(ss.nplanets+ss.nstars)),2)+' '+$
        chain_png_name + ' ' +$
        summaryfile + ' '  + summaryfile
  spawn, cmd, output


end
