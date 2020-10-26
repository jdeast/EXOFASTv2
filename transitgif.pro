;+
; NAME:
;   TRANSITGIF
;
; PURPOSE: 
;   Creates an animated gif of a transit from exoplanets.org.
;
; CALLING SEQUENCE:
;   TRANSITGIF, ['pname',PARS=,BAND=,CADENCE=,PARS=,YMIN=,BAK=,
;                /YELLOW,/DISPLAY,/SAVEFRAMES,/UPDATE]
;
; OPTIONAL INPUTS:
;   PNAME     - The name of the planet from exoplanets.org (not case
;               sensitive). The output file will be the planet name
;               stripped of spaces + .gif. 
;   PARS      - A 10 element array of custom parameters with which to
;               make a transit animation:
;                PARS[0] = period (days)
;                PARS[1] = eccentricity
;                PARS[2] = omega_star (radians)
;                PARS[3] = inclination (radians)
;                PARS[4] = a/Rstar 
;                PARS[5] = Rp/Rstar
;                PARS[6] = logg (cgs)
;                PARS[7] = Teff (K)
;                PARS[8] = [Fe/H]
;                PARS[9] = Rstar (R_sun)
;               NOTE: PNAME or PARS must be specified. If PARS is
;               specified and PNAME is not, the output will be "transit.gif".
;   BAND      - The bandpass for the limb darkening calculation (see
;               quadld.pro for allowed values). Default = 'Sloanz'
;   CADENCE   - A frame is generated every CADENCE minutes for a
;               duration of 2*T_FWHM*MAXSTAR/RSTAR. This duration
;               forces the X axis of the plot to be on the same scale as
;               the drawing, ignoring the curvature of the planet's orbit.
;               Default = 5.
;               NOTE: A constant cadence between animations means the
;               duration of the animation (number of frames)
;               corresponds to the duration of the transit.
;   NPOINTS   - The number of pixels in the image of the star. Larger
;               images may break the conversion to gif.
;               Default = 299. 
;   YMIN      - The minimum of transit plot. Default is 0.985
;               NOTE: A constant minimum makes the plots to scale from
;               animation to animation, but this default cuts off some
;               deeper planets.
;   BAK       - A 3-element array containing the RGB values of the
;               background color. Default is [135,206,250]; (sky blue)
;               Use [255,255,255] for white.
;   MAXSTAR   - The size of the maximum star, in R_sun. A star of this
;               size will be centered left/right and span the top half of
;               the image. The total image is 4*MAXSTAR
;               square. Keeping MAXSTAR the same (positive value)
;               for multiple animations means they will be to
;               scale. Default = 1.5.
;               Numbers between -1 and 0 are interpreted as a
;               fractions of RSTAR (e.g, -1 ensures the star will take
;               up the entire range; -0.9 leaves a nice margin). Negative
;               numbers are useful if you're making a single
;               animation and want it to fit nicely.
;               NOTE: If the Rstar is larger than MAXSTAR, it will
;               produce an error.
; OPTIONAL KEYWORDS:
;   YELLOW     - Makes the star yellow and the planet black, instead of
;                their color temperatures. 
;                NOTE: Output file is then PNAME.yellow.gif
;   DISPLAY    - Spawns firefox to display the animated gif.
;   SAVEFRAMES - The individual PNG frames are deleted upon
;                completion. Set this keyword to keep them.
;   UPDATE     - Update the local copy of exoplanets.csv
;
; OUTPUTS:
;  pname.gif - An animated gif of the transit.
;
; EXAMPLES:  
;   ;; creates an animation of HAT-P-3b, called 'HAT-P-3b.gif'
;   transitgif, 'HAT-P-3b'
;
;   ;; creates an animation Earth, called 'transit.gif'
;   pars = [365.242199d0,0.01671123d0,!dpi/2d0,!dpi/2d0,215.094177d0,$
;           0.0091705248d0,4.43834417d0,5778d0,0d0,1d0]
;   transitgif, pars=pars
;
;   ;; creates an animation of Jupiter, called "Jupiter.gif'
;   transitgif, 'Jupiter',pars=[4332.59d0,0.048498,!dpi/2d0,!dpi/2d0,215.094177d0*5.20260d,0.102792236d0,4.43834417d0,5778d0,0d0,1d0]
;
;   ;; Creates animations for all transiting exoplanets from exoplanets.org
;   planets = readexo(/update)
;   transit = where(planets.transit,ntransits)
;   for i=0, ntransits-1 do $
;      transitgif, strcompress(planets.name[transit[i]],/remove_all),$
;      maxstar=max(double(planets.rstar[transit])), band='Sloanz',$
;      ymin=1d0-max(double(planets.depth[transit]))
;
; DEPENDENCIES:
;   convert. While IDL can write animated gifs directly,
;   they're only 8 bit. In color, that creates unsightly
;   ringing, so we use convert instead.
;   EXOFAST (http://astroutils.astronomy.ohio-state.edu/exofast/)
;
; MODIFICATION HISTORY
;   2012/06 -- Public Release, Jason Eastman (LCOGT)
;-
pro transitgif, pname, band=band, cadence=cadence, pars=pars, ymin=ymin,$
                yellow=yellow, display=display, saveframes=saveframes,$
                update=update, bak=bak, maxstar=maxstar, npoints=npoints, $
                earth=earth, jupiter=jupiter,showpath=showpath,transparent=transparent

if keyword_set(earth) then $
   pars = [365.242199d0,0.01671123d0,!dpi/2d0,!dpi/2d0,215.094177d0,$
           0.0091705248d0,4.43834417d0,5778d0,0d0,1d0]

if keyword_set(jupiter) then $
   pars = [4332.59d0,0.048498,!dpi/2d0,!dpi/2d0,215.094177d0*5.20260d0,$
           0.102792236d0,4.43834417d0,5778d0,0d0,1d0]

if n_elements(maxstar) eq 0 then maxstar = 1.5d0
if maxstar lt -1 then message, 'MAXSTAR must be greater than -1'

if n_elements(ymin) eq 0 then ymin=0.985d0

if n_elements(band) eq 0 then band = 'Sloanz'

;; sampling rate (points/day = one every 5 minutes)
if n_elements(cadence) eq 0 then cadence = 5d0

gamma = 1440d0/cadence

if n_elements(npoints) eq 0 then npoints = 299
xsize = (npoints + 1)*2
ysize = (npoints + 1)*2

;; prepare the plotting device
mydevice = !d.name
set_plot, 'Z'
device, set_resolution=[xsize,ysize], set_pixel_depth=24
black = '000000'x

;; skyblue in RGB
if n_elements(bak) ne 3 then bak = [135,206,250]

bgcolor = 256L^2*bak[2] + 256L*bak[1] + bak[0]

;; extras not found in exoplanets.org
if n_elements(pars) eq 10 then begin
   if n_elements(pname) eq 0 then pname = 'transit'
endif else begin
   ;; get the parameters of the system from exoplanets.org
   data = readexo(update=update)
   match = (where(strupcase(strcompress(data.name,/remove_all)) eq $
                  strupcase(strcompress(    pname,/remove_all))))(0)
   if match eq -1 then begin
      print, 'pname must be one of: ' + transpose(data.name)
      return
   endif
   
   logg = double(data.logg[match])
   if logg eq 0d0 then begin
      mstar = double(data.mstar[match])
      rstar = double(data.rstar[match])
      G = 27437.4768d0
      logg = alog10(G*mstar/rstar^2)
   endif

   pars = [double(data.per[match]),double(data.ecc[match]),$
           data.om[match]*!dpi/180.d0,double(data.i[match])*!dpi/180d0,$
           double(data.ar[match]),sqrt(double(data.depth[match])),$
           logg,double(data.teff[match]),$
           double(data.fe[match]),double(data.rstar[match])]
   
   ;; check for sufficient information
   ndx = [0,3,4,5,6,7,9]
   if min(pars[ndx]) eq 0d0 then begin
      print, 'ERROR: not enough information in exoplanets.org for ' + pname
      print, pars
      return
   endif

endelse

period      = pars[0]
e           = pars[1]
omega       = pars[2]
inclination = pars[3]
ar          = pars[4]
p           = pars[5]
logg        = pars[6]
teff        = pars[7] 
feh         = pars[8]
rstar       = pars[9]

if maxstar lt 0 then maxstar = -rstar/maxstar
                                
if rstar gt maxstar then $
   message, 'ERROR: RSTAR (' + strtrim(rstar,2) + $
            ') is larger than MAXSTAR (' + strtrim(maxstar,2) + ')'

teq = teff*sqrt(1d0/(2d0*ar))
esinw = e*sin(omega)
sini = sin(inclination)
phase = exofast_getphase(e,omega,/primary)
tp = -period*phase

;; quadratic limb darkening
coeffs = quadld(logg, teff, feh, band)
u1 = coeffs[0]
u2 = coeffs[1]
if ~finite(u1) or ~finite(u2) then message, $
   'ERROR: Limb darkening tables do not cover this star ' + $
   'logg = ' + strtrim(logg,2) + ', Teff = ' + strtrim(teff,2) + $
   ', [Fe/H] = ' + strtrim(feh,2)

;; makes x scale same for plot and drawing
t14 = period/!dpi*asin(sqrt((1d0+p)^2 - p^2)/(sini*ar))*$
      sqrt(1d0-e^2)/(1d0+esinw)
t23 = period/!dpi*asin(sqrt((1d0-p)^2 - p^2)/(sini*ar))*$
      sqrt(1d0-e^2)/(1d0+esinw)
tfwhm = (t14+t23)/2d0
duration = 2d0*tfwhm/rstar*maxstar
nframes = round(duration*gamma)

;; odd # frames so bjd[nframes/2] = 0
if nframes/2d0 eq nframes/2 then nframes++ 
bjd = duration*(findgen(nframes)/(nframes-1d0) - 0.5d0)

;; define the color scaling
if keyword_set(yellow) then begin
   rflux = 255
   gflux = 255
   bflux = 0
   rplanet = 0
   gplanet = 0
   bplanet = 0
   filebase = strcompress(pname,/remove_all) + '.yellow.'
endif else begin
   ;; make a true color temperature star (based on Teff)
   readcol, getenv('EXOFAST_PATH') + 'colortable.txt', $
            temp, rtable,gtable,btable,/silent,comment='#'

   ;; steps of 100K -- don't bother interpolating; can't tell by eye anyway
   dummy = min(abs(temp-teff),match)
   rflux = rtable[match]
   gflux = gtable[match]
   bflux = btable[match]

   ;; color temperature of planet (scaled for intensity)
   dummy = min(abs(temp-teq),match)
   rplanet = round(rtable[match]*(teq/teff)^4)
   gplanet = round(gtable[match]*(teq/teff)^4)
   bplanet = round(btable[match]*(teq/teff)^4)
   
   filebase = strcompress(pname,/remove_all) + '.'

endelse

;; make the quadratically limb-darkened profile
mu = cos(asin(2d0*dindgen(npoints*rstar/maxstar/2d0)/(npoints*rstar/maxstar)))
ld = (1d0-u1*(1d0-mu)-u2*(1d0-mu)^2)

;; 1D profile in RGB (append background color to end
rprofile = [ld*rflux,dblarr(sqrt(2d0)*npoints/2d0)+bak[0]]
gprofile = [ld*gflux,dblarr(sqrt(2d0)*npoints/2d0)+bak[1]]
bprofile = [ld*bflux,dblarr(sqrt(2d0)*npoints/2d0)+bak[2]]

;; create the background images
rimage = dblarr(xsize,ysize/2) + bak[0]
gimage = dblarr(xsize,ysize/2) + bak[1]
bimage = dblarr(xsize,ysize/2) + bak[2]

;; convert 1D radial profile into circle, add it to the background
d = shift(dist(npoints+1),npoints/2,npoints/2)
rimage[xsize/4:3*xsize/4-1,*] = rprofile[d]
gimage[xsize/4:3*xsize/4-1,*] = gprofile[d]
bimage[xsize/4:3*xsize/4-1,*] = bprofile[d]
alpha = rimage*0d0 + 255B

;; get the planet's path across the star
b = exofast_getb(bjd, i=inclination,a=ar, tperiastron=tp,$
                 period=period,e=e,omega=omega,$
                 x=xplanet,y=yplanet,z=zplanet)

;; calculate the model light curve
z = sqrt(xplanet^2 + yplanet^2)
exofast_occultquad, z, u1, u2, p, flux

;; position of the planet
yplanet = -yplanet ;; looks prettier (and it doesn't matter)
xcen = (xplanet*rstar/maxstar + 1d0)*npoints/2.d0 + xsize/4
ycen = (yplanet*rstar/maxstar + 1d0)*npoints/2.d0

;; the X axis of the plot
hours = (bjd - bjd[nframes/2])*24d
xmin=min(hours,max=xmax)

;; reference grid for where the planet is
xgrid = (lonarr(ysize/2)+1)##lindgen(xsize)
ygrid = lindgen(ysize/2)##(lonarr(xsize)+1)

;; draw the planet's path over the master frame
if keyword_set(showpath) then begin
   for j=0, n_elements(xcen)-1 do begin
      pathndx = where((xgrid-xcen[j])^2+(ygrid-ycen[j])^2 le $
                   (0.5d0*xsize/299d0)^2)
      if pathndx[0] ne -1 then begin
         rimage[pathndx] = 0
         gimage[pathndx] = 0
         bimage[pathndx] = 0
      endif
   endfor
endif

for i=0, nframes-1 do begin
   
   ;; set the current frame to just the star
   rnow = rimage
   gnow = gimage
   bnow = bimage

   ;; plot the model
   device, set_font='Times', /tt_font
   plot, hours, flux, xstyle=9,ystyle=5,position=[0,0.15,1,0.45],$
         yrange=[ymin,1.002d],xrange=[xmin,xmax],color=black,$
         xtitle=TeXtoIDL('Time - T_C (hours)'),title=title,background=bgcolor, thick=xsize/299d0, $
         charsize=xsize/299d0, xthick=xsize/299d0, charthick=xsize/299d0, font=1
   ;; overplot the planet at the current frame
   plotsym, 0, /fill, color=256^2*bplanet + 256*gplanet + rplanet
   oplot, [hours[i]], [flux[i]], psym=8,symsize=xsize/299d0
   
   ;; draw the planet over the current frame
   pndx = where((xgrid-xcen[i])^2+(ygrid-ycen[i])^2 le $
                (p*npoints*rstar/maxstar/2d0)^2)
   if pndx[0] ne -1 then begin
      rnow[pndx] = rplanet
      gnow[pndx] = gplanet
      bnow[pndx] = bplanet
   endif

   ;; create png frames for each (direct gifs are 8 bit, look terrible)
   imagename = filebase + string(i,format='(i04)') + '.png'

   if keyword_set(transparent) then begin
      alpha = rimage*0d0 + 255B
      background = where(rnow eq bak[0] and gnow eq bak[1] and bnow eq bak[2])
      if background[0] ne -1 then alpha[background] = 0B
   endif
      
   ;; draw the star/planet
   tv, rnow, 0, ysize-npoints-1, 1
   tv, gnow, 0, ysize-npoints-1, 2 
   tv, bnow, 0, ysize-npoints-1, 3
;   tv, alpha
   ;oplot, xcen, ycen,/data
   write_png, imagename, tvrd(/true)

endfor

;; create the animated gif
spawn, 'convert -delay 10 ' + filebase + '????.png ' + filebase + 'gif'
if keyword_set(display) then spawn, 'firefox ' + filebase + 'gif &'

if not keyword_set(saveframes) then $
   file_delete, file_search(filebase + '????.png')

set_plot, mydevice

end
