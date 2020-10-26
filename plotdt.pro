pro plotdt, ss, psname=psname

if keyword_set(psname) then begin
   mydevice = !d.name
   set_plot, 'PS'
   aspect_ratio=1.5
;   xsize=10.5
;   ysize=xsize/aspect_ratio
;   ysize = xsize/aspect_ratio + (ss.ntran-1)*0.6
   !p.font=0
   device, filename=psname, /color, bits=24, encapsulated=0
   device, xsize=xsize,ysize=ysize
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
   if win_state[15] then wset, 15 $
   else window, 15, xsize=xsize, ysize=ysize, xpos=screen[0]/3d0, ypos=0
   red = '0000ff'x
   black = 'ffffff'x
   symsize = 0.5         
   charsize = 1.0
endelse
   ;; plotting ranges
   sigma = stddev(doptom.resid)
   dtrange = [-5d0*sigma,5d0*sigma]
   velplotmax = max(abs(doptom.vel))
   phase = (doptom.bjd-localtc) / period
   phaseplotmax = max(abs(phase))
   minvel = min(doptom.vel,max=maxvel)
   minphase = min(phase,max=maxphase)
   
   ;; determine the phase of ingress and egress to overplot on the figures
   sini = sin(inc)
   esinw = e*sin(omega)
   bp = ar*cosi*(1d0-e^2)/(1d0+esinw)
   t14 = period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sini*ar))*sqrt(1d0-e^2)/(1d0+esinw)  ;; eq 14, 16  
   t23 = period/!dpi*asin(sqrt((1d0-p)^2 - bp^2)/(sini*ar))*sqrt(1d0-e^2)/(1d0+esinw)  ;; eq 15, 16
   tau = (t14-t23)/2d0
   Tfwhm = t14-tau
   egressphase = Tfwhm*0.5d0 / period
   
   if mydevice eq 'PS' then darkblue=60 $
   else darkblue = 'ff0000'x 
   
   
   xDTtitle='Velocity (km/s)'
   yDTtitle='Orbital Phase'
   DTlabel = doptom.label

   ;; data
   multifast_dt_plotimage, -doptom.ccf2d, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                           IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, ytitle = yDTtitle, $
                           XTICKFORMAT='(A1)', title = DTlabel + ' Doppler Data', charsize=charsize, position = [0.05, 0.69, 0.90, 0.95]   
   if mydevice eq 'PS' then loadct, 39, /Silent
   oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
   if mydevice eq 'PS' then loadct, 0, /Silent
   multifast_dt_plotimage, -doptom.model, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                           IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                           XTICKFORMAT='(A1)', ytitle = modeltext, charsize=charsize, /noerase, position = [0.05, 0.43, 0.90, 0.69]
   
   ;; model
   lambdachar= '!9' + String("154B) + '!X'

   modeltext = string(lambdachar,lambda*180d0/!dpi, inc*180d0/!dpi,format='("Model (",a,"=",f0.1,", i=",f0.1,")")')
   if mydevice eq 'PS' then loadct, 39, /Silent
   oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
   if mydevice eq 'PS' then loadct, 0, /Silent
   multifast_dt_plotimage, -doptom.resid, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                           IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                           ytitle='Residuals', xtitle=xDTtitle, charsize=charsize, /noerase, position = [0.05, 0.18, 0.90, 0.43]
   
   ;; residuals
   if mydevice eq 'PS' then loadct, 39, /Silent
   oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=4, color=darkblue, linestyle=0
   oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=4, color=darkblue, linestyle=0
   if mydevice eq 'PS' then loadct, 0, /Silent
   
   ;; color bar
   dtscale = dtrange[0] + dindgen(256)/255*(dtrange[1] - dtrange[0])
   dtscale = [[dtscale],[dtscale]]
   multifast_dt_plotimage, dtscale, XRANGE=dtrange, YRANGE = [-1,1], $
                           IMGXRANGE=dtrange, IMGYRANGE = [-1,1], RANGE=dtrange, $
                           YTICKS=1, YTickformat='(A1)', xticklen=0.2, xtitle='Fractional Variation', $
                           charsize=charsize, /noerase, position = [0.05, 0.08, 0.90, 0.10]
   
endif
