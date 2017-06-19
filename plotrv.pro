pro plotrv, ss, psname=psname

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
   red = 254
   black = 0
   blue = 63
   green = 144
   symsize=0.1
endif else begin
;   set_plot, 'X'
   red = '0000ff'x
   black = 'ffffff'x
   blue = 'ff0000'x
   green = '00ff00'x
   symsize = 1d0
   charsize = 2
   device,window_state=win_state
   symsize=1
;   if win_state[20] eq 1 then wset, 20 $
;   else window, 20, retain=2
endelse

symbols = [0,8,4,0,8,4]
fills = [1,1,1,0,0,0]
colors = [black,green,blue,black,green,blue]

allmindate = !values.d_infinity
allmaxdate = -!values.d_infinity

for i=0, ss.ntel-1 do begin
   rv = *(ss.telescope[i].rvptrs)
   mindate = min(rv.bjd,max=maxdate)
   if mindate lt allmindate then allmindate = mindate
   if maxdate gt allmaxdate then allmaxdate = maxdate
endfor

nsteps = 1000
prettytime = allmindate + (allmaxdate-allmindate)*dindgen(nsteps)/(nsteps-1.d0)
allprettymodel = prettytime*0d0

nrvfit = n_elements(where(ss.planet.fitrv))
nx = ceil(sqrt(nrvfit))
ny = ceil(nrvfit/double(nx))
!p.multi = [0,nx,ny]

if not keyword_set(psname) then begin
   if win_state[20] eq 1 then wset, 20 $
   else window, 20, retain=2
endif

for i=0, ss.nplanets-1 do begin
   
   ;; pretty model without slope or gamma
   prettymodel = exofast_rv(prettytime,ss.planet[i].tp.value,$
                            ss.planet[i].period.value,0d0,ss.planet[i].K.value,$
                            ss.planet[i].e.value,ss.planet[i].omega.value,$
                            slope=0,$
                            rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value,a=ss.planet[i].ar.value,$
                            p=ss.planet[i].p.value,vsini=ss.star.vsini.value,$
                            lambda=ss.planet[i].lambda.value,$
                            u1=0d0, t0=0d0,deltarv=deltarv)
   allprettymodel += prettymodel

   ;; phase so primary is at 0.25
   prettyphase = (((prettytime - ss.planet[i].tc.value) mod $
                   ss.planet[i].period.value)/ss.planet[i].period.value + $
                  1.25d0) mod 1
   sorted = sort(prettyphase)
      
   allminrv = !values.d_infinity
   allmaxrv = -!values.d_infinity
   for j=0, ss.ntel-1 do begin
      rv = *(ss.telescope[j].rvptrs)
      modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value,$
                           ss.planet[i].period.value,0d0,$
                           ss.planet[i].K.value,ss.planet[i].e.value,$
                           ss.planet[i].omega.value,slope=0,$
                           rossiter=ss.planet[i].rossiter, i=ss.planet[i].i.value,a=ss.planet[i].ar.value,$
                           p=ss.planet[i].p.value,vsini=ss.star.vsini.value,$
                           lambda=ss.planet[i].lambda.value,$
                           u1=0d0, t0=0d0,deltarv=deltarv)
      minrv = min(rv.residuals-rv.err+modelrv)
      maxrv = max(rv.residuals+rv.err+modelrv)
      if minrv lt allminrv then allminrv = minrv
      if maxrv gt allmaxrv then allmaxrv = maxrv
   endfor

;   xtitle=TeXtoIDL('Phase + (T_P - T_C)/P + 0.25')
   xtitle1='!3Phase + (T_P - T_C)/P + 0.25'
   plot, [0], [0], xrange=[0,1], yrange=[allminrv,allmaxrv], $
         xtitle=xtitle1,ytitle='!3RV (m/s)', charsize=charsize, title=ss.planet[i].label
   for j=0, ss.ntel-1 do begin 
      rv = *(ss.telescope[j].rvptrs)
      modelrv = exofast_rv(rv.bjd,ss.planet[i].tp.value,$
                           ss.planet[i].period.value,0d0,$
                           ss.planet[i].K.value,ss.planet[i].e.value,$
                           ss.planet[i].omega.value,slope=0)
      time=(((rv.bjd-ss.planet[i].tc.value) mod ss.planet[i].period.value)/$
            ss.planet[i].period.value+1.25d0) mod 1
      plotsym, symbols[j], symsize, fill=fills[j], color=colors[j]
      oploterr, time, rv.residuals+modelrv, rv.err, 8

   endfor
   oplot, prettyphase[sorted], prettymodel[sorted], color=red
endfor

!p.multi=0

;; now plot all planets, unphased (slope and gamma subtracted)
if not keyword_set(psname) then begin
   if win_state[21] eq 1 then wset, 21 $
   else window, 21, retain=2
endif else begin
   device, /close
endelse
set_plot, mydevice

if 0 then begin
roundto = 10L^(strlen(strtrim(floor(allmaxdate-allmindate),2))+1L) 
bjd0 = floor(allmindate/roundto)*roundto

xtitle2='!3BJD_TDB - ' + string(bjd0,format='(i7)') 
plot, [0], [0], xrange=[allmindate,allmaxdate]-bjd0, $
      yrange=[allminrv,allmaxrv], xtitle=xtitle2,$
      ytitle='!3RV (m/s)'

;TeXtoIDL('BJD_{TDB} - ') + string(bjd0,format='(i7)'), 
oplot, prettytime-bjd0, allprettymodel, color=red

for j=0, ss.ntel-1 do begin 
   rv = *(ss.telescope[j].rvptrs)
   plotsym, symbols[j], symsize, fill=fills[j], color=colors[j]
   oploterr, rv.bjd-bjd0, rv.rv-ss.telescope[j].gamma.value, rv.err, 8
endfor
endif

;stop

end
