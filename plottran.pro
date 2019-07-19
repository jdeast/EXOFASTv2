pro plottran, ss, psname=psname

au = ss.constants.au/ss.constants.rsun ;; AU in rsun (~215)

aspect_ratio=1.5
mydevice=!d.name
if keyword_set(psname) then begin
   set_plot, 'PS'
   aspect_ratio=1.5
   xsize=10.5
   ysize=xsize/aspect_ratio
   ysize = xsize/aspect_ratio + (ss.ntran-1)*0.6
;   ysize=(xsize/aspect_ratio + (ss.ntran)*0.2) < screen[1]
   !p.font=0

   defsysv, '!GDL', exists=runninggdl  
   if runninggdl then psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.1.ps' $
   else psname0 = psname

   device, filename=psname0, /color, bits=24, encapsulated=0
   device, xsize=xsize,ysize=ysize
   loadct, 39, /silent
   red = 254
   symsize = 0.33
   black = 0
   charsize = 0.75
endif else begin
   !p.multi=0
   screen = GET_SCREEN_SIZE()
   device,window_state=win_state
   xsize = 600
   ysize=(xsize/aspect_ratio + (ss.ntran)*150) < screen[1]
   if win_state[30] then wset, 30 $
   else window, 30, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2
   red = '0000ff'x
   black = 'ffffff'x
   symsize = 0.5         
   charsize = 1.0
endelse
plotsym, 0, /fill, color=black

;; breaks for grazing transits...
depth = max(ss.planet.p.value^2)

noise = 0d0
for j=0, ss.ntran-1 do if stddev((*(ss.transit[j].transitptrs)).residuals) gt noise then noise = stddev((*(ss.transit[j].transitptrs)).residuals)
if keyword_set(noresiduals) then spacing = (depth+noise)*3d0 $
else spacing = (depth+noise)*4d0

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

;;secondary eclipse time
phase = exofast_getphase(e,omega,/primary)
phase2 = exofast_getphase(e,omega,/secondary)
ts = ss.planet.tc.value - ss.planet.period.value*(phase-phase2)

xmax = max(t14)*36d0
xmin = -xmax
xrange = [xmin,xmax]

j=0
trandata = (*(ss.transit[j].transitptrs)) 
time = (trandata.bjd - ss.planet[ss.transit[j].pndx].tc.value - ss.transit[j].epoch*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0

xmin = min(time,max=xmax)
xrange=[xmin,xmax]

maxnoise = stddev((*(ss.transit[0].transitptrs)).residuals)
minnoise = stddev((*(ss.transit[ss.ntran-1].transitptrs)).residuals)
if ss.ntran eq 1 then begin
   ymin = 1-3*maxnoise-depth
   ymax = 1+3*maxnoise
endif else begin
   ymin = 1d0 - depth - 3*minnoise - spacing/2d0
   ymax = 1d0 + 3*maxnoise + spacing*(ss.ntran - 0.5)
endelse
yrange = [ymin,ymax]

i=0
sini = sin(acos(ss.planet[i].cosi.value))
esinw = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
bp = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+esinw)
t14 = (ss.planet[i].period.value/!dpi*asin(sqrt((1d0+ss.planet[i].p.value)^2 - bp^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+esinw))*24d0

trandata = (*(ss.transit[0].transitptrs)) 
minbjd = min(trandata.bjd,max=maxbjd)
minbjd -= 0.25 & maxbjd += 0.25d0
t0 = floor(minbjd)

if maxbjd - minbjd lt 1 then begin
   xrange=[-t14,t14]
   xtitle='!3Time - T!DC!N (hrs)'   
endif else begin
   xrange = [minbjd,maxbjd]-t0
   xtitle='!3' + exofast_textoidl('BJD_{TDB}') + ' - ' + strtrim(t0,2)
endelse

;; output separate models and pretty models files for each planet
if n_elements(psname) eq 0 then begin
   base = 'tmpbase'
endif else base = file_dirname(psname) + path_sep() + file_basename(psname,'.transit.ps')

files = file_search(base + '.detrendedmodel.transit_*.planet_*.txt',count=nfiles)
if nfiles gt 0 then file_delete, files

;; position keyword required for proper bounding box
plot, [0],[0],yrange=yrange, xrange=xrange,/xstyle,/ystyle,$;position=[0.15, 0.05, 0.93, 0.93],$
      ytitle='!3Normalized flux + Constant',xtitle=xtitle

;; make a plot for each input file
for j=0, ss.ntran-1 do begin

   trandata = (*(ss.transit[j].transitptrs)) 
   band = ss.band[ss.transit[j].bandndx]
   
   minbjd = min(trandata.bjd,max=maxbjd)
   minbjd -= 0.25d0 & maxbjd += 0.25d0
   npretty = ceil((maxbjd-minbjd)*1440d0/5d0) ;; 1 per 5 minutes
   npoints = n_elements(trandata.bjd)

   ninterp = ss.transit[j].ninterp
   if ninterp gt 1 then begin
      prettytmptime = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
      prettytime = prettytmptime#(dblarr(ninterp)+1d0) + $
                   ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                    1440d0*ss.transit[j].exptime)##(dblarr(npretty)+1d0)
      prettyflux = dblarr(npretty,ninterp) + 1d0

      transitbjd = trandata.bjd#(dblarr(ninterp)+1d0) + $
                   ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                    1440d0*ss.transit[j].exptime)##(dblarr(npoints)+1d0)
      modelflux = dblarr(npoints,ninterp) + 1d0
   endif else begin
      prettytime = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
      prettyflux = dblarr(npretty) + 1d0
      transitbjd = trandata.bjd
      modelflux = dblarr(npoints) + 1d0
   endelse

;   ;; get the motion of the star due to the planet (pretty cadence)
;   junk = exofast_getb2(prettytime,inc=ss.planet.i.value,a=ss.planet.ar.value,$
;                        tperiastron=ss.planet.tp.value,$
;                        period=ss.planet.period.value,$
;                        e=ss.planet.e.value,omega=ss.planet.omega.value,$
;                        q=ss.star.mstar.value/ss.planet.mpsun.value,$
;                        x1=x1pretty,y1=y1pretty,z1=z1pretty)
;
;   ;; get the motion of the star due to the planets (data cadence)
;   junk = exofast_getb2(transitbjd,inc=ss.planet.i.value,a=ss.planet.ar.value,$
;                        tperiastron=ss.planet.tp.value,$
;                        period=ss.planet.period.value,$
;                        e=ss.planet.e.value,omega=ss.planet.omega.value,$
;                        q=ss.star.mstar.value/ss.planet.mpsun.value,$
;                        x1=x1,y1=y1,z1=z1)

   for i=0, ss.nplanets-1 do begin

      if ss.planet[i].fittran then begin


         ;; calculate the model for this planet at a high, regular cadence
         prettytmpflux = (exofast_tran(prettytime, $
                                       ss.planet[i].i.value + ss.transit[j].tiv.value, $
                                       ss.planet[i].ar.value, $
                                       ss.planet[i].tp.value + ss.transit[j].ttv.value, $
                                       ss.planet[i].period.value, $
                                       ss.planet[i].e.value,$
                                       ss.planet[i].omega.value,$
                                       ss.planet[i].p.value + ss.transit[j].tdeltav.value,$
                                       band.u1.value, $
                                       band.u2.value, $
                                       1d0, $
                                       q=ss.star.mstar.value/ss.planet[i].mpsun.value, $
                                       thermal=band.thermal.value, $
                                       reflect=band.reflect.value, $
                                       dilute=band.dilute.value,$
                                       tc=ss.planet[i].tc.value,$
                                       rstar=ss.star.rstar.value/AU,$
                                       ;x1=x1pretty,y1=y1pretty,z1=z1pretty,$
                                       au=au,$
                                       c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
         prettytmpflux = reform(prettytmpflux,npretty,ninterp)
         prettyflux += prettytmpflux

         ;; calculate the model for this planet for each data point
         tmpflux = (exofast_tran(transitbjd, $
                                 ss.planet[i].i.value + ss.transit[j].tiv.value, $
                                 ss.planet[i].ar.value, $
                                 ss.planet[i].tp.value + ss.transit[j].ttv.value, $
                                 ss.planet[i].period.value, $
                                 ss.planet[i].e.value,$
                                 ss.planet[i].omega.value,$
                                 ss.planet[i].p.value + ss.transit[j].tdeltav.value,$
                                 band.u1.value, $
                                 band.u2.value, $
                                 1d0, $
                                 q=ss.star.mstar.value/ss.planet[i].mpsun.value, $
                                 thermal=band.thermal.value, $
                                 reflect=band.reflect.value, $
                                 dilute=band.dilute.value,$
                                 tc=ss.planet[i].tc.value,$
                                 rstar=ss.star.rstar.value/AU,$
                                 ;x1=x1,y1=y1,z1=z1,$
                                 au=au,$
                                 c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
         tmpflux = reform(tmpflux,n_elements(trandata.bjd),ninterp)
         modelflux += tmpflux

         if ninterp gt 1 then begin
            exofast_forprint, prettytmptime, total(prettytmpflux,2)/ninterp, format='(f0.10,x,f0.10)', textout=base + '.prettymodel.transit_' + strtrim(j,2) + '.planet_' + strtrim(i,2) + '.txt', /nocomment,/silent
            exofast_forprint, trandata.bjd, total(tmpflux,2)/ninterp, format='(f0.10,x,f0.10)', textout=base + '.detrendedmodel.transit_' + strtrim(j,2) + '.planet_' + strtrim(i,2) + '.txt', /nocomment,/silent
         endif else begin
            exofast_forprint, prettytime, prettytmpflux, format='(f0.10,x,f0.10)', textout=base + '.prettymodel.transit_' + strtrim(j,2) + '.planet_' + strtrim(i,2) + '.txt', /nocomment,/silent
            exofast_forprint, transitbjd, tmpflux, format='(f0.10,x,f0.10)', textout=base + '.detrendedmodel.transit_' + strtrim(j,2) + '.planet_' + strtrim(i,2) + '.txt', /nocomment,/silent
         endelse
         
         minepoch = floor((minbjd-ss.planet[i].tc.value)/ss.planet[i].period.value)
         maxepoch = ceil((maxbjd-ss.planet[i].tc.value)/ss.planet[i].period.value)
         epochs = -minepoch + dindgen(maxepoch-minepoch+1)
         tcs = ss.planet[i].tc.value + epochs*ss.planet[i].period.value
         ;xyouts, tcs-t0, epochs*0d0+(ymax+1d0)/2d0, ss.planet[i].label, align=0.5d0
         
      endif
   endfor

   ;; now integrate the model points (before detrending)
   if ninterp gt 1 then begin
      modelflux = total(modelflux,2)/ninterp
      prettyflux = total(prettyflux,2)/ninterp
   endif

   period = ss.planet[ss.transit[j].pndx].period.value

   if maxbjd - minbjd lt period then begin
      time = (trandata.bjd - ss.planet[ss.transit[j].pndx].tc.value - ss.transit[j].epoch[ss.transit[j].pndx]*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0

      good = where(time ge xrange[0] and time le xrange[1], ngood)
      if ngood eq 0 then begin
         time = (trandata.bjd - ts[ss.transit[j].pndx] - ss.transit[j].epoch[ss.transit[j].pndx]*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0
         prettytime = (prettytime - ts[ss.transit[j].pndx] - ss.transit[j].epoch[ss.transit[j].pndx]*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0         
      endif else prettytime = (prettytime - ss.planet[ss.transit[j].pndx].tc.value  - ss.transit[j].epoch[ss.transit[j].pndx]*ss.planet[ss.transit[j].pndx].period.value + ss.transit[j].ttv.value)*24.d0

   endif else begin
      time = trandata.bjd - t0
      prettytime -= t0
   endelse

   oplot, time, modelflux + trandata.residuals + spacing*(ss.ntran-j-1), psym=8, symsize=symsize
   oplot, prettytime, prettyflux + spacing*(ss.ntran-j-1), thick=2, color=red;, linestyle=0
   xyouts, 0, 1d0 + 2d0*noise + spacing*(ss.ntran-j-1), trandata.label,charsize=charsize,alignment=0.5

endfor

;; make a phased plot for each planet
ntranfit = n_elements(where(ss.planet.fittran))
nx = ceil(sqrt(ntranfit))
ny = ceil(ntranfit/double(nx))
!p.multi = [0,nx,ny]
ysize = xsize/aspect_ratio
if keyword_set(psname) then begin
   if runninggdl then begin
      device, /close
      psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.2.ps'
      device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
   endif else device, xsize=xsize,ysize=ysize
endif else begin
   if win_state[31] then wset, 31 $
   else window, 31, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2
endelse
trandata = (*(ss.transit[0].transitptrs)) 

for i=0L, ss.nplanets-1 do begin
   
   if ss.planet[i].fittran then begin
      ymax = 1d0 + 3*maxnoise ;+ spacing*((ss.ntran-1d0)>0)

      ;; read in the detrended models from above
      prettyfiles = file_search(base + '.prettymodel.transit_*.planet_' + strtrim(i,2) + '.txt',count=nfiles)
      files = file_search(base + '.detrendedmodel.transit_*.planet_' + strtrim(i,2) + '.txt',count=nfiles)
      for j=0, nfiles-1 do begin
         readcol, files[j], thistime, thisflux, format='d,d', /silent
         readcol, prettyfiles[j], thisprettytime, thisprettyflux, format='d,d', /silent
         if j eq 0 then begin
            time = thistime-ss.transit[j].ttv.value
            modelflux = thisflux+1d0
            residuals = (*(ss.transit[j].transitptrs)).residuals
            prettytime = thisprettytime-ss.transit[j].ttv.value
            prettymodelflux = thisprettyflux+1d0
         endif else begin
            time = [time,thistime-ss.transit[j].ttv.value]
            modelflux = [modelflux,thisflux+1d0]
            residuals = [residuals,(*(ss.transit[j].transitptrs)).residuals]
            prettytime = [prettytime,thisprettytime-ss.transit[j].ttv.value]
            prettymodelflux = [prettymodelflux,thisprettyflux+1d0]
         endelse
      endfor

      ymin = min(modelflux) - 3d0*minnoise
      
      phasetime = ((time - ss.planet[i].tc.value) mod ss.planet[i].period.value)*24d0
      toohigh = where(phasetime gt (ss.planet[i].period.value/2d0*24d0))
      if toohigh[0] ne -1 then phasetime[toohigh] -= ss.planet[i].period.value*24d0
      toolow = where(phasetime lt (-ss.planet[i].period.value/2d0*24d0))
      if toolow[0] ne -1 then phasetime[toolow] += ss.planet[i].period.value*24d0
      sorted = sort(phasetime)

      prettyphasetime = ((prettytime - ss.planet[i].tc.value) mod ss.planet[i].period.value)*24d0
      toohigh = where(prettyphasetime gt (ss.planet[i].period.value/2d0*24d0))
      if toohigh[0] ne -1 then prettyphasetime[toohigh] -= ss.planet[i].period.value*24d0
      toolow = where(prettyphasetime lt (-ss.planet[i].period.value/2d0*24d0))
      if toolow[0] ne -1 then prettyphasetime[toolow] += ss.planet[i].period.value*24d0
      prettysorted = sort(prettyphasetime)

      sini = sin(acos(ss.planet[i].cosi.value))
      esinw = ss.planet[i].e.value*sin(ss.planet[i].omega.value)
      bp = ss.planet[i].ar.value*ss.planet[i].cosi.value*(1d0-ss.planet[i].e.value^2)/(1d0+esinw)
      t14 = ((ss.planet[i].period.value/!dpi*asin(sqrt((1d0+ss.planet[i].p.value)^2 - bp^2)/(sini*ss.planet[i].ar.value))*sqrt(1d0-ss.planet[i].e.value^2)/(1d0+esinw))*24d0) > (29.425d0/60d0)

      ;; plot the shell, phased model, and phased data
      ;; this plot may have some wiggles in it because of the
      ;; different ld parameters between transits
      plot, [0],[0], xstyle=1,ystyle=1,$
            ytitle='!3Normalized flux',yrange=[ymin,ymax],xrange=[-t14,t14],$
            xtitle='!3Time - Tc (Hrs)';,title=ss.planet[i].label
      oplot, phasetime, residuals + modelflux, psym=8, symsize=symsize
      oplot, prettyphasetime[prettysorted], prettymodelflux[prettysorted], thick=2, color=red, linestyle=0

   endif
endfor
!p.multi=0

;; clean up files
if n_elements(psname) eq 0 then $
   file_delete, file_search(base + '.detrendedmodel.transit_*.planet_*.txt'), /allow_nonexistent


if keyword_set(psname) then begin
   device, /close
endif
set_plot, mydevice

end
