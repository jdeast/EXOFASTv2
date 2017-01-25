;; Given a mass, [Fe/H], Teff, rstar, and uteff, this will return the
;; chi^2 penalty between the input Teff and YY Teff, and age

function massradius_yy, mstar, feh, teff, rstar, uteff,afe=afe,age=age,debug=debug,yyteff=yyteff,psname=psname, ufeh=ufeh

if n_elements(afe) eq 0 then afe = 0d0

;; convert [Fe/H] to Z (Table 2, Yi et al, 2001)
z = [[0.00001d0,-3.29d0],$ 
     [0.00010d0,-2.29d0],$  
     [0.00040d0,-1.69d0],$  
     [0.00100d0,-1.29d0],$  
     [0.00400d0,-0.68d0],$ 
     [0.00700d0,-0.43d0],$  
     [0.01000d0,-0.27d0],$  
     [0.02000d0, 0.05d0],$  
     [0.04000d0, 0.39d0],$  
     [0.06000d0, 0.60d0],$  
     [0.08000d0, 0.78d0]]
yyz = interpol(z[0,*],z[1,*],feh)

;; interpolate the YY tracks
track = yytrack(mstar, yyz, afe)

;; exclude tracks older than the universe
youngenough = where(track[0,*] le 15,nages)
if nages eq 0 then return, !values.d_infinity 
track = track[*,youngenough]

;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
sigmab = 5.670373d-5/3.839d33*6.9566d10^2 
yyteffall = 10d0^track[1,*]
yyrstarall = sqrt(10d0^track[2,*]/(4d0*!dpi*yyteffall^4d0*sigmaB))
yyageall = track[0,*]

;; locate all monotonic segments (inspired by lclxtrem.pro)
deriv = yyteffall[1:nages-1] - yyteffall[0:nages-2]
pos = where(deriv gt 0d0,complement=neg)
if pos[0] ne -1 then deriv[pos] = 1
zero = where(deriv[neg] eq 0d0)
if neg[0] ne -1 then deriv[neg] = -1
if zero[0] ne -1 then deriv[neg[zero]] = 0
deriv2 = deriv[1:nages-2] - deriv[0:nages-3]
ndx = where(deriv2 ne 0) + 1
if n_elements(ndx) eq 1 then ndx = [0,nages-1] $
else ndx = [0,ndx,nages-1]

;; if the track never crosses through this rstar, return infinity
yyteff = !values.d_infinity 
mindiff = !values.d_infinity
if arg_present(age) then age = !values.d_infinity

for i=0, n_elements(ndx)-2 do begin

   ;; do not extrapolate, only interpolate
   if ((yyrstarall[ndx[i]] ge rstar) and (yyrstarall[ndx[i+1]] le rstar)) or $
      ((yyrstarall[ndx[i]] le rstar) and (yyrstarall[ndx[i+1]] ge rstar)) then begin

      yyteffnew = interpol(yyteffall[ndx[i]:ndx[i+1]], yyrstarall[ndx[i]:ndx[i+1]],rstar)
      ;; if this segment is provides a better match to the observed Teff, use it.
      if abs(yyteffnew - teff) lt mindiff then begin
         yyteff = yyteffnew
         mindiff = abs(yyteff - teff)
         
         if arg_present(age) then $
            age = interpol(yyageall[ndx[i]:ndx[i+1]],yyrstarall[ndx[i]:ndx[i+1]],rstar)
         
      endif
   endif
endfor

chi2 = ((yyteff-teff)/uteff)^2

;; prepare the plotting device
if keyword_set(psname) then begin
   ;; astrobetter.com tip on making pretty IDL plots
   mydevice=!d.name
   set_plot, 'PS'
   aspect_ratio=1
   xsize=10.5/1.5
   ysize=xsize/aspect_ratio
   !p.font=0
   device, filename=psname, /color, bits=24
   device, xsize=xsize,ysize=ysize
   loadct, 39, /silent
   red = 254
   symsize = 0.33
endif else begin
   if keyword_set(debug) then begin
      red = '0000ff'x
      symsize = 1
      device,window_state=win_state
      if win_state[0] eq 1 then wset, 0 $
      else window, 0, retain=2
   endif
endelse

if keyword_set(debug) or keyword_set(psname) then begin

   npoints = 100
   mstar2 = 0.4 + 5d0*dindgen(npoints)/(npoints-1d0)
   yylogg2 = dblarr(npoints)
   yyteff2 = dblarr(npoints)

track = yytrack(mstar, yyz, afe)
yyteffall = 10d0^track[1,*]
yyrstarall = sqrt(10d0^track[2,*]/(4d0*!dpi*yyteffall^4d0*sigmaB))
yyageall = track[0,*]

logg = alog10(27443.4141d0*mstar/rstar^2)

   for i=0, npoints-1 do begin
      track = yytrack(mstar2[i],yyz,afe)

      yyteffall = 10d0^track[1,*]
      yyrstarall = sqrt(10d0^track[2,*]/(4d0*!dpi*yyteffall^4d0*sigmaB))
      yyageall = track[0,*]
      yyloggall = alog10(27443.4141d0*mstar2[i]/yyrstarall^2)

      dummy = min(abs(track[0,*]-age),match)
      yylogg2[i] = yyloggall[match]
      yyteff2[i] = yyteffall[match]

   endfor

   xmin=min(alog10([6000,3000]),max=xmax)
   yrange=[3,5]
;   xrange=[3.8,3.5]
   plotsym,0,/fill
   plot, alog10(yyteff2), yylogg2,xtitle=textoidl('log(T_{eff})'),ytitle='log g', xrange=[xmax,xmin],yrange=yrange
   oplot, alog10([teff]), [logg], psym=8,symsize=0.5
   oploterror, alog10([teff]), [logg], 0.00222,0.01,/lobar
   oploterror, alog10([teff]), [logg], 0.00316978,0.011,/hibar

if 0 then begin

   for j=0, 1 do begin
;      if j eq 0 then yyz2 = interpol(z[0,*],z[1,*],feh-ufeh) $
;      else yyz2 = interpol(z[0,*],z[1,*],feh+ufeh)
      
      for i=0, npoints-1 do begin
         track = yytrack(mstar2[i],yyz2,afe)
         
         yyteffall = 10d0^track[1,*]
         yyrstarall = sqrt(10d0^track[2,*]/(4d0*!dpi*yyteffall^4d0*sigmaB))
         yyageall = track[0,*]
         yyloggall = alog10(27443.4141d0*mstar2[i]/yyrstarall^2)
         
         dummy = min(abs(track[0,*]-age),match)
         yylogg2[i] = yyloggall[match]
         yyteff2[i] = yyteffall[match] 

      endfor
;      oplot, alog10(yyteff2), yylogg2,linestyle=1
   endfor
endif
;stop
if 0 then begin

   print, mstar, rstar, feh, teff, yyteff, uteff, chi2
   yyloggall = alog10(27443.4141d0*mstar/yyrstarall^2)
   logg = alog10(27443.4141d0*mstar/rstar^2)

   plotsym,0,/fill
   xmin=min(alog10(yyteffall),max=xmax)
   plot, alog10(yyteffall), yyloggall,xtitle=textoidl('log(T_{eff})'),ytitle='log g', xrange=[xmax,xmin]
   oplot, alog10([teff]), [logg], psym=8
stop
forprint, yyteffall, yyloggall, yyageall,/t
endif

;   oplot, [yyteff],[rstar],psym=2
;   oplot, yyteffall[minndx],yyrstarall[minndx],psym=2,color='00ff00'x
;   oplot, yyteffall[maxndx],yyrstarall[maxndx],psym=2,color='0000ff'x
endif

if keyword_set(psname) then begin
   device, /close
   !p.font=-1
   set_plot, mydevice
endif


return, chi2



stop




end
