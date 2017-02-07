pro plottran, ss, psname=psname
AU = 215.094177d0

aspect_ratio=1.5
if keyword_set(psname) then begin
   mydevice = !d.name
   set_plot, 'PS'
   xsize = 10.5
   ysize=xsize/aspect_ratio + (ss.ntran)*2
   !p.font=0
   device, xsize=xsize,ysize=ysize,/color,bits=24
   device, filename=psname
   loadct, 39, /silent
   red = 254
   symsize = 0.33
   black = 0
endif else begin
   !p.multi=0
   screen = GET_SCREEN_SIZE()
   device,window_state=win_state
   xsize = 600
   ysize=(xsize/aspect_ratio + (ss.ntran)*150) < screen[1]
   if win_state(3) then wset, 3 $
   else window, 3, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
   red = '0000ff'x
   black = 'ffffff'x
   symsize = 0.5         
   charsize = 1.0

endelse
plotsym, 0, /fill, color=black


;; breaks for grazing transits...
depth = max(ss.planet.p.value^2)
depth = 0.0004

if keyword_set(noresiduals) then spacing = depth*3d0 $;0.035d0 $
else spacing = depth*4d0;0.045d0



;xrange = [-3.25,3.25]
;yrange = [0.98d0,1.02d0+spacing*((ss.ntran-1d0)>0)]

cosi = ss.planet.cosi.value
sini = sin(acos(cosi))
e = ss.planet.secosw.value^2 + ss.planet.sesinw.value^2
circular = where(e eq 0d0, complement=eccentric)
omega = e*0d0
if circular[0] ne -1 then omega[circular] = !dpi/2d0 
if eccentric[0] ne -1 then omega[eccentric] = atan((ss.planet.sesinw.value)[eccentric],(ss.planet.secosw.value)[eccentric])

esinw = e*sin(omega)
ecosw = e*cos(omega)
ar = ss.planet.ar.value
period = ss.planet.period.value
p = ss.planet.p.value

bp = ar*cosi*(1d0-e^2)/(1d0+esinw)
t14 = period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sini*ar))*$
      sqrt(1d0-e^2)/(1d0+esinw)

xmax = max(t14)*36d0
xmin = -xmax
xrange = [xmin,xmax]


j=0
t0 = 2457000
trandata = (*(ss.transit[j].transitptrs)) 
time = (trandata.bjd - ss.planet[ss.transit[j].pndx].tc.value - ss.transit[j].epoch*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0
time = trandata.bjd - t0

xmin = min(time,max=xmax)
xrange=[xmin,xmax]


noise = 0.0001d0
ymin = 1d0 - 3d0*max(depth)
ymax = 1d0 + 3*noise + spacing*((ss.ntran-1d0)>0)
yrange = [ymin,ymax]
;yrange = [0.994,1.001]
xtitle='Time - T!DC!N (hrs)'
xtitle='BJD_TDB - ' + strtrim(t0,2)

;; position keyword required for proper bounding box
plot, [0],[0],yrange=yrange, xrange=xrange,/xstyle,/ystyle,position=[0.15, 0.05, 0.93, 0.93],$
      ytitle='Normalized flux + Constant',xtitle=xtitle

xrange = [xmin,xmax]

for j=0, ss.ntran-1 do begin

   trandata = (*(ss.transit[j].transitptrs)) 
   band = ss.band[ss.transit[j].bandndx]
   
   minbjd = min(trandata.bjd,max=maxbjd)
   minbjd -= 0.25 & maxbjd += 0.25d0
   npretty = ceil((maxbjd-minbjd)*1440d0) ;; 1 per minute
   prettytime = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
   prettyflux = dblarr(npretty) + 1

   for i=0, ss.nplanets-1 do begin

      if ss.planet[i].fittran then begin

         prettyflux += (exofast_tran(prettytime, $
                                     ss.planet[i].i.value, $
                                     ss.planet[i].ar.value, $
                                     ss.planet[i].tp.value, $
                                     ss.planet[i].period.value, $
                                     ss.planet[i].e.value,$
                                     ss.planet[i].omega.value,$
                                     ss.planet[i].p.value,$
                                     band.u1.value, $
                                     band.u2.value, $
                                     1d0, $
                                     thermal=band.thermal.value, $
                                     reflect=band.reflect.value, $
                                     dilute=band.dilute.value,$
                                     tc=ss.planet[i].tc.value,$
                                     rstar=ss.star.rstar.value/AU)-1d0)
         
         minepoch = floor((ss.planet[i].tc.value - minbjd)/ss.planet[i].period.value)
         maxepoch = ceil((maxbjd - ss.planet[i].tc.value)/ss.planet[i].period.value)
         epochs = -minepoch + dindgen(maxepoch+minepoch)
         tcs = ss.planet[i].tc.value + epochs*ss.planet[i].period.value
         xyouts, tcs-t0, epochs*0d0+(ymax+1d0)/2d0, ss.planet[i].label, align=0.5d0

      endif
   endfor
   prettyflux *= ss.transit[j].f0.value
   period = ss.planet[ss.transit[j].pndx].period.value

;   dummy = min(abs(trandata.bjd - ss.planet[*].tc.value),match)
   time = (trandata.bjd - ss.planet[ss.transit[j].pndx].tc.value - ss.transit[j].epoch*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0
   time = trandata.bjd - 2457000d0

;   if max(time) gt period*12d0 or min(time) lt -period*12d0 then time = time mod (period*24d0)


;   prettytime = (prettytime - ss.planet[ss.transit[j].pndx].tc.value  - ss.transit[j].epoch*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0
   prettytime = prettytime - 2457000d0

;   if max(prettytime) gt period*12d0 or min(prettytime) lt -period*12d0 then prettytime = prettytime mod (period*24d0)
;   sorted = sort(prettytime)
;   prettytime = prettytime[sorted]
;   prettyflux = prettyflux[sorted]

   oplot, time, trandata.flux + spacing*(ss.ntran-j-1), psym=8, symsize=symsize
   oplot, prettytime, prettyflux + spacing*(ss.ntran-j-1), thick=2, color=red;, linestyle=0
   xyouts, 0, 1.0075 + spacing*(ss.ntran-j-1), trandata.label,charsize=charsize,alignment=0.5
;         stop

;   wset, 1
;   plot, time, trandata.flux-
endfor

if keyword_set(psname) then begin
   device, /close
   set_plot, mydevice
endif

;for j=0, ss.nplanet-1 do begin
;endfor


end
