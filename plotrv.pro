pro plotrv, ss, psname=psname, ndx=ndx

if n_elements(ndx) eq 0 then begin
   if ss.nsteps eq 1 then ndx = 0 $
   else minchi2 = min(*ss.chi2,ndx)
endif

mydevice=!d.name
if keyword_set(psname) then begin
   set_plot, 'PS'
   aspect_ratio=1.5
   xsize=10.5
   ysize=xsize/aspect_ratio
   !p.font=0
   device, filename=psname, /color, bits=24
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

for i=0, ss.ntel-1 do begin
   rv = *(ss.telescope[i].rvptrs)
   mindate = min(rv.bjd,max=maxdate)
   if mindate lt allmindate then allmindate = mindate
   if maxdate gt allmaxdate then allmaxdate = maxdate
endfor
roundto = 10L^(strlen(strtrim(floor(allmaxdate-allmindate),2)));+1L) 
bjd0 = floor(allmindate/roundto)*roundto

allmindate = (allmindate-bjd0)*0.95 + bjd0
allmaxdate = (allmaxdate-bjd0)*1.05 + bjd0

nsteps = 1000
prettytime = allmindate + (allmaxdate-allmindate)*dindgen(nsteps)/(nsteps-1.d0)
allprettymodel = prettytime*0d0


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
   t0 = (maxtime+mintime)/2.d0
   rv.residuals = rv.rv - (ss.telescope[j].gamma.value[ndx] + (rv.bjd-t0)*ss.star.slope.value[ndx])

   for i=0, ss.nplanets-1 do begin

      if ss.planet[i].fitrv then begin
         modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value[ndx],$
                              ss.planet[i].period.value[ndx],0d0,$
                              ss.planet[i].K.value[ndx],ss.planet[i].e.value[ndx],$
                              ss.planet[i].omega.value[ndx],slope=0,$
                              rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                              p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                              lambda=ss.planet[i].lambda.value[ndx],$
                              u1=0d0, t0=0d0,deltarv=deltarv)
         ;; re-populate the residual array
         rv.residuals -= modelrv
         *(ss.telescope[j].rvptrs) = rv
      endif
   endfor
endfor

for i=0, ss.nplanets-1 do begin
   
   if ~ss.planet[i].fitrv then continue

   ;; pretty model without slope or gamma
   prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value[ndx],$
                            ss.planet[i].period.value[ndx],0d0,ss.planet[i].K.value[ndx],$
                            ss.planet[i].e.value[ndx],ss.planet[i].omega.value[ndx],$
                            slope=0,$
                            rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value[ndx],a=ss.planet[i].ar.value[ndx],$
                            p=abs(ss.planet[i].p.value[ndx]),vsini=ss.star.vsini.value[ndx],$
                            lambda=ss.planet[i].lambda.value[ndx],$
                            u1=0d0, t0=0d0,deltarv=deltarv)
   allprettymodel += prettymodel

   ;; phase so primary is at 0.25
   prettyphase = (((prettytime - ss.planet[i].tc.value[ndx]) mod $
                   ss.planet[i].period.value[ndx])/ss.planet[i].period.value[ndx] + $
                  1.25d0) mod 1
   sorted = sort(prettyphase)
      
   allminrv = min(allprettymodel,max=allmaxrv)
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
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
      t0 = (maxtime+mintime)/2.d0

      minrv = min(rv.residuals-err+modelrv)
      maxrv = max(rv.residuals+err+modelrv)

      if minrv lt allminrv then allminrv = minrv
      if maxrv gt allmaxrv then allmaxrv = maxrv
   endfor

   xtitle1='!3' + exofast_textoidl('Phase + (T_P - T_C)/P + 0.25',font=font)
   plot, [0], [0], xrange=[0,1], yrange=[allminrv,allmaxrv], $
         ytitle='!3RV (m/s)', charsize=charsize, title=title[i],position=position1,xtickformat='(A1)'
   for j=0, ss.ntel-1 do begin 
      rv = *(ss.telescope[j].rvptrs)
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
   if ss.ntel gt 1 then exofast_legend, ss.telescope.label, color=colors[indgen(ss.ntel) mod ncolors],/top,/right,psym=symbols[indgen(ss.ntel) mod nsymbols], /useplotsym, charsize=0.5, fill=fills[indgen(ss.ntel) mod nfills]

   ;; plot the residuals below
   ymin = !values.d_infinity
   ymax = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
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
      err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
      time=(((rv.bjd-ss.planet[i].tc.value[ndx]) mod ss.planet[i].period.value[ndx])/$
            ss.planet[i].period.value[ndx]+1.25d0) mod 1
      plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
      oploterr, time, rv.residuals, err, 8
   endfor
   
endfor

!p.multi=0

;; now plot all planets, unphased (slope and gamma subtracted)
if not keyword_set(psname) then begin
   if win_state[21] eq 1 then wset, 21 $
   else window, 21, retain=2
endif

xtitle2='!3' + exofast_textoidl('BJD_{TDB} - ' + string(bjd0,format='(i7)'),font=font)
plot, [0], [0], xrange=[allmindate,allmaxdate]-bjd0,/xstyle,$
      yrange=[allminrv,allmaxrv], $
      ytitle='!3RV (m/s)', position=position1,xtickformat='(A1)'

oplot, prettytime-bjd0, allprettymodel, color=red

for j=0, ss.ntel-1 do begin 
   rv = *(ss.telescope[j].rvptrs)
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
   oploterr, rv.bjd-bjd0, rv.rv-ss.telescope[j].gamma.value[ndx], err, 8
endfor
if ss.ntel gt 1 then exofast_legend, ss.telescope.label, color=colors[indgen(ss.ntel) mod ncolors],/top,/right,psym=symbols[indgen(ss.ntel) mod nsymbols], /useplotsym, charsize=0.5, fill=fills[indgen(ss.ntel) mod nfills]

;; plot the residuals below
plot, [0],[0], position=position2, /noerase, $
      xrange=[allmindate,allmaxdate]-bjd0, xtitle=xtitle2,$
      yrange=[ymin,ymax], ytitle='O-C (m/s)', $
      /xstyle, /ystyle, yminor=2,yticks=2, ytickv=[ymin*0.7,0,ymax*0.7]
oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

for j=0L, ss.ntel-1 do begin
   rv = *(ss.telescope[j].rvptrs)
   err = sqrt(rv.err^2 + ss.telescope[j].jittervar.value[ndx])
   plotsym, symbols[j mod nsymbols], symsize, fill=fills[j mod nfills], color=colors[j mod ncolors]
   oploterr, rv.bjd-bjd0, rv.residuals, err, 8
endfor

if keyword_set(psname) then begin
   device, /close
   set_plot, mydevice
endif

end
