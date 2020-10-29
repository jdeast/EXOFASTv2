pro plottran, ss, psname=psname, ndx=ndx, noresiduals=noresiduals

if n_elements(ndx) eq 0 then ndx = 0L

au = ss.constants.au/ss.constants.rsun ;; AU in rsun (~215)

aspect_ratio=1.5
mydevice=!d.name
if keyword_set(psname) then begin
   set_plot, 'PS'
   aspect_ratio=1.5
   xsize=10.5
   ysize=xsize/aspect_ratio
   ysize = xsize/aspect_ratio; + (ss.ntran-1)*0.6
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
   position1 = [0.23, 0.40, 0.95, 0.95]    ;; data plot
   position2 = [0.23, 0.20, 0.95, 0.40]    ;; residual plot

endif else begin
   !p.multi=0
   !p.font=-1
   screen = GET_SCREEN_SIZE()
   device,window_state=win_state
   if win_state[25] then wset, 25 $
   else window, 25, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2
   xsize = 600
   ysize=(xsize/aspect_ratio + (ss.ntran)*150) < screen[1]
   red = '0000ff'x
   black = 'ffffff'x
   symsize = 0.5         
   charsize = 1.0
   position1 = [0.07, 0.22, 0.97, 0.95]    ;; data plot
   position2 = [0.07, 0.07, 0.97, 0.22]    ;; residual plot

endelse
plotsym, 0, /fill, color=black

;; derive some quantities to scale the plot with
tc = ss.planet.tc.value[ndx]
cosi = ss.planet.cosi.value[ndx]
sini = sin(acos(cosi))
e = ss.planet.e.value[ndx]
omega = ss.planet.omega.value[ndx]
esinw = e*sin(omega)
ecosw = e*cos(omega)
ar = ss.planet.ar.value[ndx]
period = ss.planet.period.value[ndx]
p = ss.planet.p.value[ndx]
bp = ar*cosi*(1d0-e^2)/(1d0+esinw)

;; primary duration
t14 = period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sini*ar))*$
      sqrt(1d0-e^2)/(1d0+esinw)

;;secondary eclipse time, duration
phase = exofast_getphase(e,omega,/primary)
phase2 = exofast_getphase(e,omega,/secondary)
ts = ss.planet.tc.value[ndx] - ss.planet.period.value[ndx]*(phase-phase2)
bs = ar*cosi*(1d0-e^2)/(1d0-esinw)
t14s = period/!dpi*asin(sqrt((1d0+p)^2 - bs^2)/(sini*ar))*$
      sqrt(1d0-e^2)/(1d0-esinw)

;; depth for each band in each transit
noise = dblarr(ss.ntran)
bandnames = strarr(ss.ntran)
depth = dblarr(ss.nplanets, ss.ntran)
depth2 = dblarr(ss.nplanets, ss.ntran)
primary = bytarr(ss.nplanets, ss.ntran)
secondary = bytarr(ss.nplanets, ss.ntran)
phasecurve = bytarr(ss.nplanets, ss.ntran)
minbjd = !values.d_infinity
maxbjd = -!values.d_infinity
for i=0L, ss.nplanets-1 do begin
   for j=0L, ss.ntran-1 do begin
      band = ss.band[ss.transit[j].bandndx] 
      bandnames[j] = band.name
      depth2[i,j] = band.thermal.value[ndx]/1d6

      if band.reflect.fit then phasecurve[i,j] = 1B
      if band.ellipsoidal.fit then phasecurve[i,j] = 1B
      if ss.planet[i].beam.fit ne 0 then phasecurve[i,j] = 1B

      u1 = band.u1.value[ndx]
      u2 = band.u2.value[ndx]
      exofast_occultquad_cel, abs(ss.planet[i].b.value[ndx]), u1, u2, ss.planet[i].p.value[ndx],mu1
      depth[i,j] = 1d0 - mu1

      ;; noise for each band
      noise[j] = stddev((*(ss.transit[j].transitptrs)).residuals)

      ;; time
      thisbjd = (*(ss.transit[j].transitptrs)).bjd
      thisminbjd = min(thisbjd,max=thismaxbjd)
      if thisminbjd lt minbjd then minbjd = thisminbjd
      if thismaxbjd gt maxbjd then maxbjd = thismaxbjd

      if ss.planet[i].fittran then begin
         ;; does this LC have data corresponding to the primary or
         ;; secondary of this planet (within time +/- duration)?
         minepoch = floor((thisminbjd - tc[i])/period[i])-2 ;; pad to account for rounding errors
         maxepoch = ceil((thismaxbjd - tc[i])/period[i])+2  ;; pad to account for rounding errors
         epochs = minepoch + lindgen(maxepoch-minepoch+1)
         for k=0L, n_elements(epochs)-1 do begin
            match1 = where(thisbjd gt (tc[i] + epochs[k]*period[i] - t14[i]) and $
                           thisbjd lt (tc[i] + epochs[k]*period[i] + t14[i]))
            if match1[0] ne -1 then primary[i,j] = 1B
            match2 = where(thisbjd gt (ts[i] + epochs[k]*period[i] - t14s[i]) and $
                           thisbjd lt (ts[i] + epochs[k]*period[i] + t14s[i]))
            if match2[0] ne -1 then secondary[i,j] = 1B
         endfor
      endif

   endfor
endfor


maxnoise = max(noise)
roundto = 10L^(strlen(strtrim(floor(maxbjd-minbjd),2)))
t0 = floor(minbjd/roundto)*roundto
xrange = [minbjd,maxbjd]-t0
xtitle='!3' + exofast_textoidl('BJD_{TDB}') + ' - ' + strtrim(t0,2)

;; for the first plot
yrange = [1d0-max(depth)-3*maxnoise,1d0+3*maxnoise]

;; output separate models and pretty models files for each planet
if n_elements(psname) eq 0 then begin
   base = 'tmpbase'
endif else base = file_dirname(psname) + path_sep() + file_basename(psname,'.transit.ps')

files = file_search(base + '.detrendedmodel.transit_*.planet_*.txt',count=nfiles)
if nfiles gt 0 then file_delete, files

;; draw the shell of the unphased plot
;; position keyword required for proper bounding box
plot, [0],[0],yrange=yrange, xrange=xrange,/xstyle,/ystyle,$;position=[0.15, 0.05, 0.93, 0.93],$
      ytitle='!3Normalized flux',xtitle=xtitle

;; make one plot for all input files
for j=0, ss.ntran-1 do begin

   trandata = (*(ss.transit[j].transitptrs)) 
   band = ss.band[ss.transit[j].bandndx]
   
   minbjd = min(trandata.bjd,max=maxbjd)
   minbjd -= 1d0 & maxbjd += 1d0
   
   npoints = n_elements(trandata.bjd)
   if ss.debug then npretty = npoints $
   else npretty = ceil((maxbjd-minbjd)*1440d0/5d0) ;; 1 per 5 minutes

   ninterp = ss.transit[j].ninterp
   if ninterp gt 1 then begin
      transitbjd = trandata.bjd#(dblarr(ninterp)+1d0) + $
                   ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                    1440d0*ss.transit[j].exptime)##(dblarr(npoints)+1d0)
      modelflux = dblarr(npoints,ninterp) + 1d0

      prettytmptime = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
      prettytime = prettytmptime#(dblarr(ninterp)+1d0) + $
                   ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                    1440d0*ss.transit[j].exptime)##(dblarr(npretty)+1d0)
      prettyflux = dblarr(npretty,ninterp) + 1d0

   endif else begin
      transitbjd = trandata.bjd
      modelflux = dblarr(npoints) + 1d0
      if ss.debug then begin
         prettytime = transitbjd
         prettyflux = modelflux
      endif else begin         
         prettytime = minbjd + (maxbjd-minbjd)*dindgen(npretty)/(npretty-1d0)
         prettyflux = dblarr(npretty) + 1d0
      endelse
   endelse

   for i=0, ss.nplanets-1 do begin

      if ~ss.planet[i].fittran then continue

      ;; calculate the model for this planet at a high, regular cadence
      prettytmpflux = (exofast_tran(prettytime, $
                                    ss.planet[i].i.value[ndx] + ss.transit[j].tiv.value[ndx], $
                                    ss.planet[i].ar.value[ndx], $
                                    ss.planet[i].tp.value[ndx] + ss.transit[j].ttv.value[ndx], $
                                    ss.planet[i].period.value[ndx], $
                                    ss.planet[i].e.value[ndx],$
                                    ss.planet[i].omega.value[ndx],$
                                    ss.planet[i].p.value[ndx] + ss.transit[j].tdeltav.value[ndx],$
                                    band.u1.value[ndx], $
                                    band.u2.value[ndx], $
                                    1d0, $
                                    q=ss.star.mstar.value[ndx]/ss.planet[i].mpsun.value[ndx], $
                                    thermal=band.thermal.value[ndx], $
                                    reflect=band.reflect.value[ndx], $
                                    dilute=band.dilute.value[ndx],$
                                    ellipsoidal=band.ellipsoidal.value[ndx],$
                                    beam=ss.planet[i].beam.value[ndx],$
                                    tc=ss.planet[i].tc.value[ndx],$
                                    rstar=ss.star.rstar.value[ndx]/AU,$
                                    au=au,$
                                    c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
      prettytmpflux = reform(prettytmpflux,npretty,ninterp)
      prettyflux += prettytmpflux
      
      ;; calculate the model for this planet for each data point
      tmpflux = (exofast_tran(transitbjd, $
                              ss.planet[i].i.value[ndx] + ss.transit[j].tiv.value[ndx], $
                              ss.planet[i].ar.value[ndx], $
                              ss.planet[i].tp.value[ndx] + ss.transit[j].ttv.value[ndx], $
                              ss.planet[i].period.value[ndx], $
                              ss.planet[i].e.value[ndx],$
                              ss.planet[i].omega.value[ndx],$
                              ss.planet[i].p.value[ndx] + ss.transit[j].tdeltav.value[ndx],$
                              band.u1.value[ndx], $
                              band.u2.value[ndx], $
                              1d0, $
                              q=ss.star.mstar.value[ndx]/ss.planet[i].mpsun.value[ndx], $
                              thermal=band.thermal.value[ndx], $
                              reflect=band.reflect.value[ndx], $
                              dilute=band.dilute.value[ndx],$
                              ellipsoidal=band.ellipsoidal.value[ndx],$
                              beam=ss.planet[i].beam.value[ndx],$
                              tc=ss.planet[i].tc.value[ndx],$
                              rstar=ss.star.rstar.value[ndx]/AU,$
                              au=au,$
                              c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
      tmpflux = reform(tmpflux,n_elements(trandata.bjd),ninterp)
      modelflux += tmpflux
      
      planetno = string(i,format='(i02)')
      transitno = string(j,format='(i03)')
      if ninterp gt 1 then begin
         exofast_forprint, prettytmptime, total(prettytmpflux,2)/ninterp, format='(f0.10,x,f0.10)', textout=base + '.prettymodel.transit_' + transitno + '.planet_' + planetno + '.txt', /nocomment,/silent
         exofast_forprint, trandata.bjd, total(tmpflux,2)/ninterp, format='(f0.10,x,f0.10)', textout=base + '.detrendedmodel.transit_' + transitno + '.planet_' + planetno + '.txt', /nocomment,/silent
      endif else begin
         exofast_forprint, prettytime, prettytmpflux, format='(f0.10,x,f0.10)', textout=base + '.prettymodel.transit_' + transitno + '.planet_' + planetno + '.txt', /nocomment,/silent
         exofast_forprint, transitbjd, tmpflux, format='(f0.10,x,f0.10)', textout=base + '.detrendedmodel.transit_' + transitno + '.planet_' + planetno + '.txt', /nocomment,/silent
      endelse
      
      minepoch = floor((minbjd-ss.planet[i].tc.value[ndx])/ss.planet[i].period.value[ndx])
      maxepoch = ceil((maxbjd-ss.planet[i].tc.value[ndx])/ss.planet[i].period.value[ndx])
      epochs = -minepoch + dindgen(maxepoch-minepoch+1)
      tcs = ss.planet[i].tc.value[ndx] + epochs*ss.planet[i].period.value[ndx]
      
   endfor

   ;; now integrate the model points (before detrending)
   if ninterp gt 1 then begin
      modelflux = total(modelflux,2)/ninterp
      prettyflux = total(prettyflux,2)/ninterp
   endif

   time = trandata.bjd - t0
   prettytime -= t0
   oplot, time, modelflux + trandata.residuals, psym=8, symsize=symsize
   oplot, prettytime, prettyflux, thick=2, color=red;, linestyle=0

endfor

;   if win_state[30] then wset, 30 $
;   else window, 30, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2


;; make a plot for each planet, with each file phased and on top of one another
;; page 2 - 2+nplanets-1
pageno = 1

for jj=0L, 1 do begin
   
   pageno++

   if jj eq 0L then begin
      mapping = primary 
      duration = t14
      delta = depth
      xtitle = exofast_textoidl('Time - T_C (Hrs)')
      t_eclipse = ss.planet.tc.value[ndx]
   endif else begin
      mapping = secondary
      duration = t14s
      delta = depth2
      xtitle = exofast_textoidl('Time - T_S (Hrs)') 
      t_eclipse = ss.planet.ts.value[ndx]
   endelse

   ;; get the right windows/page numbers
   if keyword_set(psname) then begin
      if runninggdl then begin
         device, /close
         psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(pageno,2) + '.ps'
         device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
      endif else device, xsize=xsize,ysize=ysize
   endif else begin
      if win_state[24+pageno] then wset, 24+pageno $
      else window, 24+pageno, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2
   endelse

   for i=0L, ss.nplanets-1 do begin
      
      use = where(mapping[i,*],nlc)
      if nlc eq 0 then continue

      ysize = xsize/aspect_ratio + (nlc-1)
;      ysize = xsize/aspect_ratio + (ss.ntran-1)*0.6

      
      ;; this is a planet specific quantity
      spacing = (max(delta[i,use])+max(noise[use]))*3d0
      if ~keyword_set(noresiduals) then spacing *= 4d0/3d0
      
      xrange=[-duration[i],duration[i]]*24d0
      if nlc eq 1 then begin
         ymin = 1-3*max(noise[use])-max(delta[i,use])
         ymax = 1+3*max(noise[use])
      endif else begin
         ymin = 1d0 - max(delta[i,use]) - 3d0*noise[use[0]] - spacing/2d0
         ymax = 1d0 + 3d0*noise[use[nlc-1]] + spacing*(nlc - 0.5)
      endelse
      yrange = [ymin,ymax]
      
      ;; draw the shell of the plot
      if keyword_set(psname) then begin
         if runninggdl then begin
            device, /close
            psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(pageno,2) + '.ps'
            device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
            pageno++
         endif else device, xsize=xsize,ysize=ysize
      endif
      
      plot, [0],[0],yrange=yrange, xrange=xrange,/xstyle,/ystyle,$
            ytitle='!3Normalized flux + Constant',xtitle=xtitle
      
      npretty = ceil((2d0*duration[i])*1440d0*2d0) ;; 1 per 30 seconds
      
      ;; make a plot for each input file
      for j=0, ss.ntran-1 do begin
         
         if ~mapping[i,j] then continue
         
         trandata = (*(ss.transit[j].transitptrs)) 
         band = ss.band[ss.transit[j].bandndx]
         npoints = n_elements(trandata.bjd)
         
         ninterp = ss.transit[j].ninterp
         if ninterp gt 1 then begin
            transitbjd = trandata.bjd#(dblarr(ninterp)+1d0) + $
                         ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                          1440d0*ss.transit[j].exptime)##(dblarr(npoints)+1d0)
            modelflux = dblarr(npoints,ninterp) + 1d0
            prettytmptime = t_eclipse[i] - duration[i] + (2d0*duration[i])*dindgen(npretty)/(npretty-1d0)
            
            prettytime = prettytmptime#(dblarr(ninterp)+1d0) + $
                         ((dindgen(ninterp)/ninterp-(ninterp-1d0)/(2d0*ninterp))/$
                          1440d0*ss.transit[j].exptime)##(dblarr(npretty)+1d0)
            prettyflux = dblarr(npretty,ninterp) + 1d0
            
         endif else begin
            transitbjd = trandata.bjd
            modelflux = dblarr(npoints) + 1d0
            prettytime = t_eclipse[i] - duration[i] + (2d0*duration[i])*dindgen(npretty)/(npretty-1d0)
            prettyflux = dblarr(npretty) + 1d0
         endelse
         
         ;; calculate the model for this planet at a high, regular cadence
         prettytmpflux = (exofast_tran(prettytime, $
                                       ss.planet[i].i.value[ndx] + ss.transit[j].tiv.value[ndx], $
                                       ss.planet[i].ar.value[ndx], $
                                       ss.planet[i].tp.value[ndx] + ss.transit[j].ttv.value[ndx], $
                                       ss.planet[i].period.value[ndx], $
                                       ss.planet[i].e.value[ndx],$
                                       ss.planet[i].omega.value[ndx],$
                                       ss.planet[i].p.value[ndx] + ss.transit[j].tdeltav.value[ndx],$
                                       band.u1.value[ndx], $
                                       band.u2.value[ndx], $
                                       1d0, $
                                       q=ss.star.mstar.value[ndx]/ss.planet[i].mpsun.value[ndx], $
                                       thermal=band.thermal.value[ndx], $
                                       reflect=band.reflect.value[ndx], $
                                       dilute=band.dilute.value[ndx],$
                                       ellipsoidal=band.ellipsoidal.value[ndx],$
                                       beam=ss.planet[i].beam.value[ndx],$                                      
                                       tc=ss.planet[i].tc.value[ndx],$
                                       rstar=ss.star.rstar.value[ndx]/AU,$
                                       au=au,$
                                       c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
         prettytmpflux = reform(prettytmpflux,npretty,ninterp)
         prettyflux += prettytmpflux
         
         ;; calculate the model for this planet for each data point
         tmpflux = (exofast_tran(transitbjd, $
                                 ss.planet[i].i.value[ndx] + ss.transit[j].tiv.value[ndx], $
                                 ss.planet[i].ar.value[ndx], $
                                 ss.planet[i].tp.value[ndx] + ss.transit[j].ttv.value[ndx], $
                                 ss.planet[i].period.value[ndx], $
                                 ss.planet[i].e.value[ndx],$
                                 ss.planet[i].omega.value[ndx],$
                                 ss.planet[i].p.value[ndx] + ss.transit[j].tdeltav.value[ndx],$
                                 band.u1.value[ndx], $
                                 band.u2.value[ndx], $
                                 1d0, $
                                 q=ss.star.mstar.value[ndx]/ss.planet[i].mpsun.value[ndx], $
                                 thermal=band.thermal.value[ndx], $
                                 reflect=band.reflect.value[ndx], $
                                 dilute=band.dilute.value[ndx],$
                                 ellipsoidal=band.ellipsoidal.value[ndx],$
                                 beam=ss.planet[i].beam.value[ndx],$
                                 tc=ss.planet[i].tc.value[ndx],$
                                 rstar=ss.star.rstar.value[ndx]/AU,$
                                 au=au,$
                                 c=ss.constants.c/ss.constants.au*ss.constants.day) - 1d0) 
         tmpflux = reform(tmpflux,n_elements(trandata.bjd),ninterp)
         modelflux += tmpflux
         
         ;; why am i doing this?
         minepoch = floor((minbjd-ss.planet[i].tc.value[ndx])/ss.planet[i].period.value[ndx])
         maxepoch = ceil((maxbjd-ss.planet[i].tc.value[ndx])/ss.planet[i].period.value[ndx])
         epochs = -minepoch + dindgen(maxepoch-minepoch+1)
         tcs = ss.planet[i].tc.value[ndx] + epochs*ss.planet[i].period.value[ndx]
         
         ;; now integrate the model points (before detrending)
         if ninterp gt 1 then begin
            modelflux = total(modelflux,2)/ninterp
            prettyflux = total(prettyflux,2)/ninterp
            prettytime = total(prettytime,2)/ninterp
         endif
         
         time = (trandata.bjd - (t_eclipse[i] + ss.transit[j].ttv.value[ndx])) mod ss.planet[i].period.value[ndx]
         toohigh = where(time gt ss.planet[i].period.value[ndx]/2d0)
         toolow = where(time lt -ss.planet[i].period.value[ndx]/2d0)
         if toohigh[0] ne -1 then time[toohigh] -= ss.planet[i].period.value[ndx]
         if toolow[0] ne -1 then time[toolow] += ss.planet[i].period.value[ndx]
         time *= 24d0 ;; hours to days
         
         prettytime = (prettytime - (t_eclipse[i] + ss.transit[j].ttv.value[ndx]))*24d0
         
         factor = nlc - total(mapping[i,0:j])
         oplot, time, modelflux + trandata.residuals + spacing*factor, psym=8, symsize=symsize
         oplot, prettytime, prettyflux + spacing*factor, thick=2, color=red ;, linestyle=0
         xyouts, 0, 1d0 + 2d0*max(noise) + spacing*factor, trandata.label,charsize=charsize,alignment=0.5

      endfor
      
   endfor
endfor

;; make a phased plot for each planet
;; stack the transits on top of each other
ysize = xsize/aspect_ratio
for jj=0L, 2 do begin
   
   if jj eq 0L then begin
      ;; make a phased primary transit plot
      mapping = primary 
      duration = t14
      delta = depth
      xtitle = exofast_textoidl('Time - T_C (Hrs)')
      t_eclipse = ss.planet.tc.value[ndx]
   endif else if jj eq 1 then begin
      ;; make a phased secondary eclipse plot
      mapping = secondary
      duration = t14s
      delta = depth2
      xtitle = exofast_textoidl('Time - T_S (Hrs)') 
      t_eclipse = ss.planet.ts.value[ndx]
   endif else if jj eq 2 then begin
      ;; make a phase curve plot
      mapping = phasecurve
      duration = period/2d0
      delta = depth2
      xtitle=exofast_textoidl('Time - T_S (Hrs)')
      t_eclipse = ss.planet.ts.value[ndx]
   endif
   
   for i=0L, ss.nplanets-1 do begin
      
      if ~ss.planet[i].fittran then continue
      match = where(mapping[i,*],nlc)
      if nlc eq 0 then continue

      pageno++

      sorted = sort(bandnames[match])
      sortband = (bandnames[match])[sorted]
      uniqbands = sortband[uniq(sortband)]
      for kk=0L, n_elements(uniqbands)-1 do begin
         
         ;; get the right windows/page numbers
         if keyword_set(psname) then begin
            if runninggdl then begin
               device, /close
               psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(pageno,2) + '.ps'
               device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
            endif else device, xsize=xsize,ysize=ysize
         endif else begin
            win_num = (24+pageno) mod 32
            if win_state[win_num] then wset, win_num $
            else window, win_num, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0, retain=2
         endelse
         
         use = where(mapping[i,*] and bandnames eq uniqbands[kk],nuse)
         maxnoise = max(noise[use])
         ymax = max(modelflux)+3d0*maxnoise
         
         ;; read in the detrended models from above
         planetno = string(i,format='(i02)')
         prettyfiles = file_search(base + '.prettymodel.transit_*.planet_' + planetno + '.txt',count=nfiles)
         files = file_search(base + '.detrendedmodel.transit_*.planet_' + planetno + '.txt',count=nfiles)
         for j=0, nuse-1 do begin
            readcol, files[use[j]], thistime, thisflux, format='d,d', /silent
            readcol, prettyfiles[use[j]], thisprettytime, thisprettyflux, format='d,d', /silent
            if j eq 0 then begin
               time = thistime-ss.transit[use[j]].ttv.value[ndx]
               modelflux = thisflux+1d0
               residuals = (*(ss.transit[use[j]].transitptrs)).residuals
               prettytime = thisprettytime-ss.transit[use[j]].ttv.value[ndx]
               prettymodelflux = thisprettyflux+1d0
            endif else begin
               time = [time,thistime-ss.transit[use[j]].ttv.value[ndx]]
               modelflux = [modelflux,thisflux+1d0]
               residuals = [residuals,(*(ss.transit[use[j]].transitptrs)).residuals]
               prettytime = [prettytime,thisprettytime-ss.transit[use[j]].ttv.value[ndx]]
               prettymodelflux = [prettymodelflux,thisprettyflux+1d0]
            endelse
         endfor

         phasetime = ((time - t_eclipse[i]) mod ss.planet[i].period.value[ndx])*24d0
         toohigh = where(phasetime gt (ss.planet[i].period.value[ndx]/2d0*24d0))
         if toohigh[0] ne -1 then phasetime[toohigh] -= ss.planet[i].period.value[ndx]*24d0
         toolow = where(phasetime lt (-ss.planet[i].period.value[ndx]/2d0*24d0))
         if toolow[0] ne -1 then phasetime[toolow] += ss.planet[i].period.value[ndx]*24d0
         sorted = sort(phasetime)

         if jj eq 2 then begin
            ;; ignore the primary transit for the Y scaling (to
            ;; highlight phase variations)
            use = where(phasetime le (ss.planet[i].period.value[ndx] - t14) and $
                        phasetime ge (-ss.planet[i].period.value[ndx] + t14))
         endif else begin
            use = where(abs(phasetime) lt duration[i])
         endelse

         ymin = min(modelflux[use]) - 3d0*maxnoise

         
         prettyphasetime = ((prettytime - t_eclipse[i]) mod ss.planet[i].period.value[ndx])*24d0
         toohigh = where(prettyphasetime gt (ss.planet[i].period.value[ndx]/2d0*24d0))
         if toohigh[0] ne -1 then prettyphasetime[toohigh] -= ss.planet[i].period.value[ndx]*24d0
         toolow = where(prettyphasetime lt (-ss.planet[i].period.value[ndx]/2d0*24d0))
         if toolow[0] ne -1 then prettyphasetime[toolow] += ss.planet[i].period.value[ndx]*24d0
         prettysorted = sort(prettyphasetime)
         
         if keyword_set(psname) then begin
            if runninggdl then begin
               device, /close
               psname0 = file_dirname(psname) + path_sep() + file_basename(psname,'.ps') + '.' + strtrim(pageno,2) + '.ps'
               device, filename=psname0, /color, bits=24, xsize=xsize,ysize=ysize
               pageno++
            endif else device, xsize=xsize,ysize=ysize
         endif
         
         ;; plot the shell, phased model, and phased data
         ;; this plot may have some wiggles in it because of the
         ;; different ld parameters between transits
         plotbandname = uniqbands[kk]
         if strpos(plotbandname,'Sloan') ne -1 then $
            plotbandname = (strsplit(plotbandname,'Sloan',/regex,/extract))[0] + "'"
         
         xrange = [-duration[i],duration[i]]*24d0
         plot, [0],[0], xstyle=1,ystyle=1,$
               ytitle='!3Normalized flux (' + plotbandname + ')',yrange=[ymin,ymax],xrange=xrange,$
               position=position1, xtickformat='(A1)'
         oplot, phasetime, residuals + modelflux, psym=8, symsize=symsize
         if jj eq 2 then begin
            oplot, phasetime+ss.planet[i].period.value[ndx]*24d0, residuals + modelflux, psym=8, symsize=symsize
            oplot, phasetime-ss.planet[i].period.value[ndx]*24d0, residuals + modelflux, psym=8, symsize=symsize
         endif

         oplot, prettyphasetime[prettysorted], prettymodelflux[prettysorted], thick=2, color=red, linestyle=0
         if jj eq 2 then begin
            oplot, prettyphasetime[prettysorted]+ss.planet[i].period.value[ndx]*24d0, prettymodelflux[prettysorted], thick=2, color=red, linestyle=0
            oplot, prettyphasetime[prettysorted]-ss.planet[i].period.value[ndx]*24d0, prettymodelflux[prettysorted], thick=2, color=red, linestyle=0
         endif
         
         ;; plot the residuals below
         ymin = min([0,residuals*1.1],max=ymax)
         if ymin lt -ymax then ymax = -ymin
         if ymax gt -ymin then ymin = -ymax
         yrange = [ymin,ymax] 
         plot, [0],[0], position=position2, /noerase, $
               xrange=xrange, xtitle=xtitle,$
               yrange=yrange, ytitle='O-C', $
               /xstyle, /ystyle, yminor=2,yticks=2,ytickv=[ymin,0d0,ymax]*0.7d0
         oplot, phasetime, residuals, psym=8,symsize=symsize
         if jj eq 2 then begin
            oplot, phasetime+ss.planet[i].period.value[ndx]*24d0, residuals, psym=8,symsize=symsize
            oplot, phasetime-ss.planet[i].period.value[ndx]*24d0, residuals, psym=8,symsize=symsize
         endif
         oplot, [-9d9,9d9],[0,0],linestyle=2,color=red  

      endfor
   endfor
endfor

;; clean up files
if n_elements(psname) eq 0 then $
   file_delete, file_search(base + '.detrendedmodel.transit_*.planet_*.txt'), /allow_nonexistent


if keyword_set(psname) then begin
   device, /close
endif
set_plot, mydevice

end
