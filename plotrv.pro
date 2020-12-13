pro plotrv, ss, psname=psname, ndx=ndx, yrange=yrange, savfile=savfile

if n_elements(savfile) ne 0 then begin
   restore, savfile
   ss = mcmcss
endif

;; pick the best-fit index if not specified
if n_elements(ndx) eq 0 then begin
   if ss.nsteps eq 1 then ndx = 0 $
   else minchi2 = min(*ss.chi2,ndx)
endif

defsysv, '!GDL', exists=runninggdl
mydevice=!d.name
if keyword_set(psname) then begin
   set_plot, 'PS'
   aspect_ratio=1.5
   xsize=10.5
   ysize=xsize/aspect_ratio
   !p.font=0

   if runninggdl then psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.1.ps' $
   else psname0 = psname

   device, filename=psname0, /color, bits=24
   device, xsize=xsize,ysize=ysize
   loadct, 39, /silent
   colors = [0,159,95,254,223,31,207,111,191,47]
   red = 254
   symsize=0.5
   title = strarr(ss.nplanets)
   position1 = [0.23, 0.40, 0.95, 0.95]    ;; data plot
   position2 = [0.23, 0.20, 0.95, 0.40]    ;; residual plot
endif else begin
;   set_plot, 'X'
   colors= ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x,'0080ff'x,$
            '800080'x,'00ffff'x,'ffff00'x,'80d000'x,'660000'x]
   red = '0000ff'x
   symsize = 1d0
   charsize = 1
   device,window_state=win_state
   symsize=1

   title = ss.planet.label
   position1 = [0.07, 0.22, 0.97, 0.95]    ;; data plot
   position2 = [0.07, 0.07, 0.97, 0.22]    ;; residual plot
   font=-1

;   if win_state[20] eq 1 then wset, 20 $
;   else window, 20, retain=2
endelse

symbols = [0,8,4,3,0,8,4,3]
fills = [1,1,1,1,0,0,0,0]
nsymbols = n_elements(symbols)
nfills = n_elements(fills)
ncolors = n_elements(colors)

allmindate = !values.d_infinity
allmaxdate = -!values.d_infinity

q = ss.star.mstar.value[ndx]/ss.planet.mpsun.value[ndx]

starrvs = [-1]
companionrvs = [-1]

legendndx = lonarr(ss.nplanets+1, ss.ntel)
for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   mindate = min(rv.bjd,max=maxdate)
   if mindate lt allmindate then allmindate = mindate
   if maxdate gt allmaxdate then allmaxdate = maxdate
   
   if rv.planet eq -1 then legendndx[ss.nplanets,j] = 1B $
   else legendndx[rv.planet,j] = 1B

endfor
t0 = (allmindate+allmaxdate)/2d0

roundto = 10L^(strlen(strtrim(floor(allmaxdate-allmindate),2)));+1L) 
bjd0 = floor(allmindate/roundto)*roundto

allmindate = (allmindate-bjd0)*0.95 + bjd0
allmaxdate = (allmaxdate-bjd0)*1.05 + bjd0

;; if the timespan is > 100 years, there's got to be an error in the timestamps
if (allmaxdate-allmindate) gt 36500 then begin
   message, 'WARNING: RV data span > 100 years. Make sure your input timestamps are consistent. This may cause memory issues', /continue
endif

;; sample the pretty light curve 100 times per (minimum) period 
;; for the span of the data
cadence = min(ss.planet.period.value[ndx])/100d0
nsteps = (allmaxdate-allmindate)/cadence
prettytime = allmindate + (allmaxdate-allmindate)*dindgen(nsteps)/(nsteps-1.d0)
allprettymodel = prettytime*0d0
allprettymodel2 = prettytime*0d0

if not keyword_set(psname) then begin
   if win_state[20] eq 1 then wset, 20 $
   else window, 20, retain=2
   nrvfit = n_elements(where(ss.planet.fitrv))
   ny = ceil(sqrt(nrvfit))
   nx = ceil(nrvfit/double(ny))
   !p.multi = [0,nx,ny]
endif

;; calculate the residuals
for j=0, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   mintime = min(rv.bjd,max=maxtime)

   ;; subtract gamma, slope, and quadratic terms
   rv.residuals = rv.rv - (ss.telescope[j].gamma.value[ndx] + ss.star.slope.value[ndx]*(rv.bjd-t0) + ss.star.quad.value[ndx]*(rv.bjd-t0)^2)

   for i=0, ss.nplanets-1 do begin

      if ~ss.planet[i].fitrv then continue
      if rv.planet ne -1 then continue

      if rv.planet eq i then begin
;; this needs to be debugged
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i])
         
         ;; calculate the RV model
         modelrv += exofast_rv(rvbjd,ss.planet[i].tp.value[ndx],ss.planet[i].period.value[ndx],$
                               0d0,ss.planet[i].K.value[ndx]*q[i],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi,$
                               slope=0d0)
         
      endif else begin
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i], /primary)

         modelrv = exofast_rv(rvbjd,ss.planet[i].tp.value[ndx],$
                              ss.planet[i].period.value[ndx],0d0,$
                              ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                              ss.planet[i].omega.value[ndx],slope=0,$
                              rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                              p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                              lambda=ss.planet[i].lambda.value[ndx],$
                              u1=0d0,deltarv=deltarv)
         ;; re-populate the residual array
         rv.residuals -= modelrv
         *(ss.telescope[j].rvptrs) = rv

      endelse
   endfor
endfor

for i=0, ss.nplanets-1 do begin
   
   if ~ss.planet[i].fitrv then continue
   if rv.planet ne -1 then continue

   if rv.planet eq i then begin
;; this needs to be debugged

      ;; pretty model without quad, slope, or gamma
      prettymodel2 = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                                ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx]*q[i],$
                                ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi)
      allprettymodel2 += prettymodel2


   endif else begin

      ;; pretty model without quad, slope, or gamma
      prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                               ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx],$
                               slope=0,$
                               rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                               p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                               lambda=ss.planet[i].lambda.value[ndx],$
                               u1=0d0, deltarv=deltarv)
      allprettymodel += prettymodel

      if n_elements(psname) eq 1 then begin
         base = file_dirname(psname) + path_sep() + file_basename(psname,'.model')
         exofast_forprint, prettytime, prettymodel, textout=base+'.prettymodelrv.planet.' + strtrim(i,2) + '.txt', format='(f0.10,x,f0.10)'
      endif
   endelse 

   ;; phase so primary is at 0.25
   prettyphase = (((prettytime - ss.planet[i].tc.value[ndx]) mod $
                   ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                  1.25d0) mod 1
   sorted = sort(prettyphase)
   
   if n_elements(yrange) eq 0 then begin
      allminrv = min(allprettymodel,max=allmaxrv)
      for j=0, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         if rv.planet ne -1 then continue

         err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
         modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx],$
                              ss.planet[i].period.value[ndx],0d0,$
                              ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                              ss.planet[i].omega.value[ndx],slope=0,$
                              rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                              p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                              lambda=ss.planet[i].lambda.value[ndx],$
                              u1=0d0, t0=0d0,deltarv=deltarv)
         
         mintime = min(rv.bjd,max=maxtime)
         minrv = min(rv.residuals-err+modelrv)
         maxrv = max(rv.residuals+err+modelrv)
         
         if minrv lt allminrv then allminrv = minrv
         if maxrv gt allmaxrv then allmaxrv = maxrv
      endfor
   endif else begin
      allminrv = yrange[0]
      allmaxrv = yrange[1]
   endelse
   
   xtitle1='!3' + exofast_textoidl('Phase + (T_P - T_C)/P + 0.25',font=font)
   plot, [0], [0], xrange=[0,1], yrange=[allminrv,allmaxrv], $
         ytitle='!3RV (m/s)', charsize=charsize, title=title[i],position=position1,xtickformat='(A1)'
   for j=0, ss.ntel-1 do begin 
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx],$
                           ss.planet[i].period.value[ndx],0d0,$
                           ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                           ss.planet[i].omega.value[ndx],slope=0,$
                           rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                           p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                           lambda=ss.planet[i].lambda.value[ndx],$
                           u1=0d0, t0=0d0,deltarv=deltarv)
      time=(((rv.bjd-ss.planet[i].tc.value[ndx]) mod ss.planet[i].period.value[ndx])/$
            ss.planet[i].period.value[ndx]+1.25d0) mod 1
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals+modelrv, err, 8

   endfor
   oplot, prettyphase[sorted], prettymodel[sorted], color=red
   use = where(legendndx[ss.nplanets,*],nuse)
   if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/bottom,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]

   ;; plot the residuals below
   ymin = !values.d_infinity
   ymax = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      ymin = min([(rv.residuals - err)*1.1,ymin])
      ymax = max([(rv.residuals + err)*1.1,ymax])
   endfor

   ;; make the plot symmetric about 0
   if ymin lt -ymax then ymax = -ymin
   if ymax gt -ymin then ymin = -ymax

   plot, [0],[0], position=position2, /noerase, $
         xrange=xrange, xtitle=xtitle1,$
         yrange=[ymin,ymax], ytitle='O-C (m/s)', $
         /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin*0.7,0,ymax*0.7]
   oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
   
   for j=0L, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      if rv.planet ne -1 then continue
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      time=(((rv.bjd-ss.planet[i].tc.value[ndx]) mod ss.planet[i].period.value[ndx])/$
            ss.planet[i].period.value[ndx]+1.25d0) mod 1
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals, err, 8
   endfor

   ;; GDL can't do multi-page plots
   if runninggdl then begin
      device, /close
      psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(i+2,2) + '.ps'
      device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
   endif
endfor

!p.multi=0

;; now plot all planets, unphased, including the slope and quadratic terms
if not keyword_set(psname) then begin
   if win_state[21] eq 1 then wset, 21 $
   else window, 21, retain=2
endif else begin
   trend = (prettytime-t0)*ss.star.slope.value[ndx] + (prettytime-t0)^2*ss.star.quad.value[ndx]
   allprettymodel += trend
   exofast_forprint, prettytime, prettymodel, textout=base+'.prettymodelrv.trend.txt', format='(f0.10,x,f0.10)'
endelse

if n_elements(yrange) eq 0 then begin
;; scale for the unphased plot
   allminrv = min(allprettymodel,max=allmaxrv)
   for j=0, ss.ntel-1 do begin
      rvstr = *(ss.telescope[j].rvptrs)
      if rvstr.planet ne -1 then continue
      rv = rvstr.rv - ss.telescope[j].gamma.value[ndx]
      err = sqrt(rvstr.err^2 + ss.telescope[j].jittervar.value[ndx])
      if min(rv-err) lt allminrv then allminrv = min(rv-err)
      if max(rv+err) gt allmaxrv then allmaxrv = max(rv+err)
   endfor
endif else begin
   allminrv = yrange[0]
   allmaxrv = yrange[1]
endelse

xtitle2='!3' + exofast_textoidl('BJD_{TDB} - ' + string(bjd0,format='(i7)'),font=font)
plot, [0], [0], xrange=[allmindate,allmaxdate]-bjd0,/xstyle,$
      yrange=[allminrv,allmaxrv], $
      ytitle='!3RV (m/s)', position=position1,xtickformat='(A1)'
oplot, prettytime-bjd0, allprettymodel, color=red

for j=0, ss.ntel-1 do begin 
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then continue
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
   oploterr, rv.bjd-bjd0, rv.rv-ss.telescope[j].gamma.value[ndx], err, 8
endfor

use = where(legendndx[ss.nplanets,*],nuse)
if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/top,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]

;; plot the residuals below
plot, [0],[0], position=position2, /noerase, $
      xrange=[allmindate,allmaxdate]-bjd0, xtitle=xtitle2,$
      yrange=[ymin,ymax], ytitle='O-C (m/s)', $
      /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin*0.7,0,ymax*0.7]
oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then continue
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
   oploterr, rv.bjd-bjd0, rv.residuals, err, 8
endfor

;; plot the companion RVs
nplanetrvs = 0L
for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   if rv.planet ne -1 then nplanetrvs++
endfor

if nplanetrvs gt 0L then begin

   for i=0L, ss.nplanets-1 do begin
      
      prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                               ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx]*q[i],$
                               ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi)
      
      ;; phase so primary is at 0.25
      prettyphase = (((prettytime - ss.planet[i].tc.value[ndx]) mod $
                      ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                     1.25d0) mod 1
      sorted = sort(prettyphase)
      
      
      allminphase = 1d0
      allmaxphase = 0d0
      allminrv = !values.d_infinity
      allmaxrv = -!values.d_infinity
      allminoc = !values.d_infinity
      allmaxoc = -!values.d_infinity
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         if rv.planet ne i then continue
         
         ;; this needs to be debugged
         rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                            a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                            period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                            omega=ss.planet[i].omega.value[ndx]+!dpi,$
                            c=ss.constants.c/ss.constants.au*ss.constants.day,$
                            q=q[i])
         
         phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                       ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                      1.25d0) mod 1
         
         minphase = min(phasetime,max=maxphase)
         if minphase lt allminphase then allminphase = minphase
         if maxphase gt allmaxphase then allmaxphase = maxphase
         minrv = min(rv.rv,max=maxrv)
         if minrv lt allminrv then allminrv = minrv
         if maxrv gt allmaxrv then allmaxrv = maxrv
      endfor
      
;   xtitle1='!3' + exofast_textoidl('Phase + (T_P - T_C)/P + 0.25',font=font)
;   xtitle2='!3' + exofast_textoidl('BJD_{TDB} - ' + string(bjd0,format='(i7)'),font=font)
;   plot, [0], [0], xtitle=xtitle1, psym=8, ytitle='!3RV (m/s)', xrange=xrange, yrange=yrange
      
      plot, [0], [0], xrange=[allminphase,allmaxphase],/xstyle,$
            yrange=[allminrv,allmaxrv], $
            ytitle='!3RV (m/s)', position=position1,xtickformat='(A1)'
      
      oplot, prettyphase, prettymodel,  color=red
      
      use = where(legendndx[i,*],nuse)
      if nuse gt 1 then exofast_legend, ss.telescope[use].label, color=colors[(indgen(ss.ntel) mod ncolors)[use]],/bottom,/right,psym=symbols[(indgen(ss.ntel) mod nsymbols)[use]], /useplotsym, charsize=0.5, fill=fills[(indgen(ss.ntel) mod nfills)[use]]
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         if rv.planet eq i then begin
            
            ;; this needs to be debugged
            rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                               a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                               period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                               omega=ss.planet[i].omega.value[ndx]+!dpi,$
                               c=ss.constants.c/ss.constants.au*ss.constants.day,$
                               q=q[i])
            
            ;; calculate the RV model
            modelrv = exofast_rv(rvbjd,ss.planet[i].tp.value[ndx],ss.planet[i].period.value[ndx],$
                                 0d0,ss.planet[i].K.value[ndx]*q[i],$
                                 ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx]+!dpi,$
                                 slope=0d0)
            
            ;; for residual plot
            (*(ss.telescope[j].rvptrs)).residuals = rv.rv-modelrv-ss.telescope[j].gamma.value[ndx]
            minoc = min((*(ss.telescope[j].rvptrs)).residuals ,max=maxoc)
            
            if minoc lt allminoc then allminoc = minoc
            if maxoc gt allmaxoc then allmaxoc = maxoc
            
            phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                          ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                         1.25d0) mod 1
            
            plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
            oplot, phasetime, modelrv, psym=8
            
         endif
      endfor
      
      ;; **** this is a hack; do this better ****
      if finite(allminoc) and finite(allmaxoc) then begin
         ;; residual plot
         plot, [0],[0], position=position2, /noerase, $
               xrange=[allminphase,allmaxphase], xtitle=xtitle1,$
               yrange=[allminoc,allmaxoc], ytitle='O-C (m/s)', $
               /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[allminoc*0.7,0,allmaxoc*0.7]
         oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  
      endif
      
      for j=0L, ss.ntel-1 do begin
         rv = *(ss.telescope[j].rvptrs)
         
         if rv.planet eq i then begin
            
            rvbjd = bjd2target(rv.bjd, inclination=ss.planet[i].i.value[ndx], $
                               a=ss.planet[i].a.value[ndx], tp=ss.planet[i].tp.value[ndx], $
                               period=ss.planet[i].period.value[ndx], e=ss.planet[i].e.value[ndx],$
                               omega=ss.planet[i].omega.value[ndx]+!dpi,$
                               c=ss.constants.c/ss.constants.au*ss.constants.day,$
                               q=q[i])
            phasetime = ((( rvbjd - ss.planet[i].tc.value[ndx]) mod $
                          ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                         1.25d0) mod 1
            plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
            oplot , phasetime, rv.residuals, psym=8
         endif
      endfor
      
   endfor
endif

if keyword_set(psname) then begin
   device, /close
   set_plot, mydevice
endif

end
