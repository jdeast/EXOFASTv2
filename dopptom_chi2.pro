;; Calculate the Doppler Tomography model
;; Written 2015 (?) by Thomas G. Beatty
;; Modified to integrate into EXOFASTv2 2017/06/27, JDE
;; note vsini and vline are in km/s here!

function dopptom_chi2, doptom, tc, period, e, omega, cosi, p, ar, lambda, logg, teff, feh, vsini, vline, errscale, debug=debug,like=like,psname=psname

if vline lt 0d0 then return, !values.d_infinity ; invalid line width
if vsini le 0d0 then return, !values.d_infinity ; invalid vsini

chisqr = 0.0d0;
convolve_limit = 5. ;; The gaussian broadening has to be less than this factor larger than vsini for it to try and include the rotation kernel. E.g., gaussian broadening ('GaussTerm') = 10km/s, vsini needs to be at least 2km/s. Otherwise just use a gaussian for the shadow.

inc = acos(cosi)
;b = ar * cosi

velsini = doptom.vel/vsini
stepsize = doptom.stepsize / vsini
meanstepsize = mean(stepsize)

c = 299792.458d0
rvel = c/doptom.rspec ;; spectrograph resolution in velocity 

;; convert form FWHM to sigma
fwhm2sigma = 2d0*sqrt(2d0*alog(2d0))

;relevantVels = where((velsini gt -1.5) and (velsini lt 1.5))  ; we only care about these, to speed things up
relevantVels = where(abs(velsini) le 2d0)  ; we only care about these, to speed things up
;relevantVels = lindgen(n_elements(velsini))

IndepVels = (rvel/fwhm2sigma) / (meanstepsize*vsini)     ; accounts for the fact that we're supersampling within the actual spectral res. of the spectrograph, to adjust chi^2 later

GaussTerm = sqrt(vline^2. + rvel^2.)/fwhm2sigma ; so we can add in both the broadening from the instrument's spectral resolution, and the vline parameter

GaussRel = GaussTerm/vsini

;; init for planet shadow model
velwidth = p
c1 = 2d0 / (!dpi*velwidth)

;xp = ar * sin(phase*2d0*!dpi)
;zp = ar * cos(phase*2d0*!dpi) * cosi
;up = xp*cos(lambda) + zp*sin(lambda)

sz = size(doptom.ccf2d)
nvels = sz[1]
ntimes = sz[2]

;; Calculate Collier-Cameron's "beta" as the actual transit LC in V
;; This presumes we're using an ~optical spectrograph for the observations
beta = dblarr(ntimes)
    
ldcoeffs = quadld(logg, teff, feh, 'V')
u1 = ldcoeffs[0]
u2 = ldcoeffs[1]
if ((not finite(u1)) or (not finite(u2))) then  return, !values.d_infinity

tp = tc - period*exofast_getphase(e,omega,/pri)
z = exofast_getb2(doptom.bjd, i=inc, a=ar, tperiastron=tp, period=Period, e=e,omega=omega,z=depth, x=xp,y=yp)
up = xp*cos(lambda) - yp*sin(lambda)

primary = where(depth lt 0,complement=secondary)     
if primary[0] ne -1 then begin
   exofast_occultquad_cel, z[primary], u1, u2, p, mu1
   beta[primary] = mu1
endif

beta = 1d0-beta                   ; to just get the amount of light blocked, rather than the ligthcurve

doptom.model = dblarr(nvels,ntimes) + median(doptom.ccf2d)
nrelvel = n_elements(relevantVels)

for i = 0, ntimes-1 do begin
   if beta[i] gt 0.0d0 then begin                                ; only do this if the planet is on the star
      c2 = ((velsini[relevantVels]-up[i])/velwidth)^2.
      validlocs = where(c2 lt 1.0d0)
      ;; Need to check if we have any valid locations
      if validlocs[0] ne -1 then begin
         rotprofile = dblarr(nrelvel)
         rotprofile[validlocs] = c1*sqrt(1d0-c2[validlocs])
         unnormalized = gaus_convol(velsini[relevantVels], rotprofile, GaussRel)
         normalization = 1d0 / total(unnormalized*stepsize)
         doptom.model[relevantVels,i] += (beta[i] * normalization * unnormalized)
      endif
   endif
endfor
model = doptom.model

if 0 then begin
good = where(beta gt 0d0,ngood)
c2 = ((velsini[relevantVels]#replicate(1d0,ngood) - up[good]##replicate(1d0,nrelvel))/velwidth)^2
validlocs = where(c2 lt 1d0)
G = dblarr(nvels,ntimes)
if validlocs[0] ne -1 then G[validlocs] = c1*sqrt(1d0-c2[validlocs])

unnormalized = dblarr(nrelvel,ngood)
for i=0, ngood-1 do unnormalized[*,good[i]] = gaus_convol(velsini[relevantVels],G[*,good[i]],GaussRel)
normalization = 1d0/total(unnormalized*(stepsize#replicate(1d0,ntimes)),1)
doptom.model[relevantvels,good] = ((beta*normalization)##replicate(1d0,nrelvel))*unnormalized
stop
endif

;; divide chi^2 by the number of independent velocities to account for correlated errors
doptom.resid = doptom.ccf2d - doptom.model
if n_elements(errscale) eq 1 then $
   chi2 = exofast_like(doptom.resid,0d0,doptom.rms*errscale,/chi2)/IndepVels $
else chi2 = total((doptom.resid/doptom.rms)^2)/IndepVels

if ~finite(chi2) then stop

if not keyword_set(psname) and not keyword_set(debug) then return, chi2

mydevice = !d.name
aspect_ratio=1.5
if keyword_set(psname) then begin
   set_plot, 'PS'
   !p.font=0
   device, filename=psname, /color, bits=24, encapsulated=1
   device, xsize=14,ysize=24;,xoffset=5.0,yoffset=2.0,/portrait
   loadct, 39, /silent
   red = 254
   symsize = 0.33
   black = 0
   darkblue=60 
   charsize=1.3
   thick=4
endif else begin
   !p.multi=0
   screen = GET_SCREEN_SIZE()
   device,window_state=win_state
   if win_state[15] then wset, 15 $
   else window, 15;, xpos=screen[0]/3d0, ypos=0
   red = '0000ff'x
   black = 'ffffff'x
   symsize = 0.5         
   charsize = 1.0
   darkblue = 'ff0000'x 
   thick=4
endelse

;; plotting ranges
sigma = stddev(doptom.resid)
dtrange = [-5d0*sigma,5d0*sigma]
velplotmax = max(abs(doptom.vel))

nper = round((mean(doptom.bjd)-tc) / period)
localtc = tc + Period*nper
phase = (doptom.bjd-localtc) / period

phaseplotmax = max(abs(phase))
minvel = min(doptom.vel,max=maxvel)
minphase = min(phase,max=maxphase)

;; determine the phase of ingress and egress to overplot on the figures
sini = sin(inc)
esinw = e*sin(omega)
bp = ar*cosi*(1d0-e^2)/(1d0+esinw)
t14 = period/!dpi*asin(sqrt((1d0+p)^2 - bp^2)/(sini*ar))*sqrt(1d0-e^2)/(1d0+esinw)    ;; eq 14, 16  
t23 = period/!dpi*asin(sqrt((1d0-p)^2 - bp^2)/(sini*ar))*sqrt(1d0-e^2)/(1d0+esinw)    ;; eq 15, 16
tau = (t14-t23)/2d0
Tfwhm = t14-tau
egressphase = Tfwhm*0.5d0 / period

xDTtitle='Velocity (km/s)'
yDTtitle='Orbital Phase'
DTlabel = doptom.label

x1 = 0.2
x2 = 0.9

;; data
if keyword_set(psname) then loadct, 0, /silent
multifast_dt_plotimage, -doptom.ccf2d, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                        IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, ytitle = yDTtitle, $
                        XTICKFORMAT='(A1)', title = DTlabel + ' Doppler Data', charsize=charsize, position = [x1, 0.69, x2, 0.95],/ystyle
if keyword_set(psname) then loadct, 39, /silent
oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=thick, color=darkblue, linestyle=0



;; model
lambdachar= '!9' + String("154B) + '!X'
modeltext = string(lambdachar,lambda*180d0/!dpi, inc*180d0/!dpi,format='("Model (",a,"=",f0.1,", i=",f0.1,")")')

if keyword_set(psname) then loadct, 0, /silent
multifast_dt_plotimage, -doptom.model, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                        IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                        XTICKFORMAT='(A1)', ytitle = modeltext, charsize=charsize, /noerase, position = [x1, 0.43, x2, 0.69],/ystyle
if keyword_set(psname) then loadct, 39, /silent
oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=thick, color=darkblue, linestyle=0

;; residuals
if keyword_set(psname) then loadct, 0, /silent
multifast_dt_plotimage, -doptom.resid, XRANGE = [-velplotmax,velplotmax], YRANGE = [-phaseplotmax,phaseplotmax], $
                        IMGXRANGE = [minvel,maxvel], IMGYRANGE = [minphase,maxphase], RANGE = DTrange, $
                        ytitle='Residuals', xtitle=xDTtitle, charsize=charsize, /noerase, position = [x1, 0.18, x2, 0.43],/ystyle
if keyword_set(psname) then loadct, 39, /silent
oplot, [-vsini,-vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [ vsini, vsini],[-phaseplotmax,phaseplotmax], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[-egressphase,-egressphase], thick=thick, color=darkblue, linestyle=0
oplot, [-velplotmax, velplotmax],[ egressphase, egressphase], thick=thick, color=darkblue, linestyle=0
   
;; color bar
dtscale = dtrange[0] + dindgen(256)/255*(dtrange[1] - dtrange[0])
dtscale = [[dtscale],[dtscale]]
if keyword_set(psname) then loadct, 0, /silent
multifast_dt_plotimage, dtscale, XRANGE=dtrange, YRANGE = [-1,1], $
                        IMGXRANGE=dtrange, IMGYRANGE = [-1,1], RANGE=dtrange, $
                        YTICKS=1, YTickformat='(A1)', xticklen=0.2, xtitle='Fractional Variation', $
                        charsize=charsize, /noerase, position = [x1, 0.08, x2, 0.10]
   
if keyword_set(psname) then begin
   !p.font=-1
   !p.multi=0
   device, /close
   device, encapsulated=0
endif
set_plot, mydevice

return, chi2

end
