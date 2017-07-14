function massradius_yy3, mstar, feh, age, teff, yyrstar=yyrstar,uteff,afe=afe,$
                         debug=debug,yyteff=yyteff, quad=quad, lsquad=lsquad,$
                         yytrack=yytrack, psname=psname, $
                         sigmab=sigmab, gravitysun=gravitysun

if n_elements(sigmab) eq 0 then $   
   sigmab = 5.670367d-5/3.828d33*6.957d10^2 ;; Stefan-boltzmann Constant (L_sun/(r_sun^2*K^4))
if n_elements(solargravity) eq 0 then $
   gravitysun = 27420.011d0 ;; cm/s^2

;; this is the approximate error in the YY isochrone modeling
uteff = 50d0 
if n_elements(afe) eq 0 then afe = 0d0

;; make these global for speed (restore is expensive)
COMMON YY_BLOCK, tracks

path = getenv('EXOFAST_PATH') + 'yy/'

if file_test(path + 'tracks.idl') then begin
   if n_elements(tracks) eq 0 then restore, path + 'tracks.idl' 
   sz = size(tracks)
   na = sz[1]
   nz = sz[2]
   nx = sz[3]
endif else begin
   ;; generate the tracks array and save it
   aname = ['a0o2','a2o2','a4o2']
   zname = ['x76997z00001', 'x7697z0001', 'x7688z0004',$
            'x767z001', 'x758z004', 'x749z007', 'x74z01',$
            'x71z02', 'x65z04', 'x59z06', 'x53z08']
   xmname = ['m04','m05','m06','m07','m08','m09','m10',$
             'm11','m12','m13','m14','m15','m16','m17','m18','m19','m20',$
             'm21','m22','m23','m24','m25','m26','m27','m28','m29','m30',$
             'm32','m34','m36','m38','m40','m42','m45','m50']
   na = n_elements(aname)
   nz = n_elements(zname)
   nx = n_elements(xmname)

   ;; nafe x nz x nmass x nage x 3 (teff, rstar, age) array 
   track1 = dblarr(na,nz,nx,24,3) + !values.d_nan
   track2 = dblarr(na,nz,nx,150,3) + !values.d_nan
   tracks = dblarr(na,nz,nx,174,3) + !values.d_nan


   for i=0,na-1 do begin
      for j=0,nz-1 do begin
         for k=0,nx-1 do begin 
            filename = path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j]
            
            ;; read all the first tracks
            if file_test(filename + '.track1') then begin
               readcol, path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j] + '.track1',$
                        n, age1,logteff1,loglstar1,format='i,d,d,d',skip=1,/silent
               track1[i,j,k,*,0] = 10^logteff1
               track1[i,j,k,*,1] = sqrt(10d0^loglstar1/(4d0*!dpi*track1[i,j,k,*,0]^4d0*sigmaB))
               track1[i,j,k,*,2] = age1
;               plot, logteff, loglstar
            endif
            
            ;; read all the second tracks
            if file_test(filename + '.track2') then begin
               readcol, path + aname[i] + path_sep() + zname[j] + path_sep() + xmname[k] + zname[j] + '.track2',$
                        n, age2,logteff2,loglstar2,format='i,d,d,d',skip=1,/silent
               track2[i,j,k,*,0] = 10^logteff2
               track2[i,j,k,*,1] = sqrt(10d0^loglstar2/(4d0*!dpi*track2[i,j,k,*,0]^4d0*sigmaB))
               track2[i,j,k,*,2] = age2
            endif

            if file_test(filename + '.track1') and file_test(filename + '.track2') then begin
               logteff = [logteff1,logteff2]
               loglstar = [loglstar1,loglstar2]
               age = [age1,age2]
               sorted = sort(age)
               tracks[i,j,k,*,0] = 10^logteff[sorted]
               tracks[i,j,k,*,1] = sqrt(10d0^loglstar[sorted]/(4d0*!dpi*tracks[i,j,k,*,0]^4d0*sigmaB))
               tracks[i,j,k,*,2] = age[sorted]
;               plot, tracks[i,j,k,*,0], tracks[i,j,k,*,1]
;               stop
            endif


         endfor
      endfor
   endfor
   save, filename=path + 'tracks.idl',tracks
endelse

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

yyz = interpol(z[0,*],z[1,*],feh,quad=quad, lsquad=lsquad)

avalue = [0d0,0.3d0,0.6d0]
zvalue = [0.00001d0,0.0001d0,0.0004d0,0.001d0,0.004d0,$
          0.007d0,0.01d0,0.02d0,0.04d0,0.06d0,0.08d0]
xmass = [0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,$
         1.1d0,1.2d0,1.3d0,1.4d0,1.5d0,1.6d0,1.7d0,1.8d0,1.9d0,2.0d0,$
         2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,2.6d0,2.7d0,2.8d0,2.9d0,3.0d0,$
         3.2d0,3.4d0,3.6d0,3.8d0,4.0d0,4.2d0,4.5d0,5.0d0]

z_one = interpol(indgen(nz),zvalue,yyz,quad=quad, lsquad=lsquad)
a_one = interpol(indgen(na),avalue,afe,quad=quad, lsquad=lsquad)
m_one = interpol(indgen(nx),xmass,mstar,quad=quad, lsquad=lsquad)

;; if the track never crosses through this rstar, return infinity
yyteff = !values.d_infinity 
mindiff = !values.d_infinity
;age = !values.d_infinity

;if keyword_set(debug) then print, mstar, rstar

;; out of range
if z_one gt (nz-1) or a_one gt (na-1) or m_one gt (nx-1) or $
   z_one lt 0 or a_one lt 0 or m_one lt 0 then return, !values.d_infinity

;; this saves a lot of time (and imposes a solid prior) by not
;; interpolating/allowing values older than the universe
yytrack = dblarr(3,sz[4])
nages = 0

if afe eq 0d0 then begin
   repeat begin
      ;; 2D cubic convolution interpolation (metalicity and mass)
      yytrack[0,nages] = interpolate(tracks[0,*,*,nages,0],z_one,m_one,cubic=-0.5d0) ;; teff
      yytrack[1,nages] = interpolate(tracks[0,*,*,nages,1],z_one,m_one,cubic=-0.5d0) ;; rstar
      yytrack[2,nages] = interpolate(tracks[0,*,*,nages,2],z_one,m_one,cubic=-0.5d0) ;; age
      
      if ~finite(yytrack[0,nages]) then return, !values.d_infinity
      if ~finite(yytrack[1,nages]) then return, !values.d_infinity
      if ~finite(yytrack[2,nages]) then return, !values.d_infinity
      
      nages++
   endrep until yytrack[2,nages-1] ge 13.9 or nages eq sz[4]
endif else begin
   repeat begin
      ;; 3D trilinear interpolation (afe, metalicity, and mass)
      yytrack[0,nages] = interpolate(tracks[*,*,*,nages,0],a_one,z_one,m_one) ;; teff
      yytrack[1,nages] = interpolate(tracks[*,*,*,nages,1],a_one,z_one,m_one) ;; rstar
      yytrack[2,nages] = interpolate(tracks[*,*,*,nages,2],a_one,z_one,m_one) ;; age

      if ~finite(yytrack[0,nages]) then return, !values.d_infinity
      if ~finite(yytrack[1,nages]) then return, !values.d_infinity
      if ~finite(yytrack[2,nages]) then return, !values.d_infinity

      nages++
   endrep until yytrack[2,nages-1] ge 15 or nages eq sz[4]
endelse

yytrack = yytrack[*,0:nages-1]

;; interpolate at the given age
yyteff  = interpol(yytrack[0,*], yytrack[2,*],age, quad=quad, lsquad=lsquad)
yyrstar = interpol(yytrack[1,*], yytrack[2,*],age, quad=quad, lsquad=lsquad)

if yyrstar le 0 then return, !values.d_infinity
;chi2 = ((yyteff-teff)/uteff)^2
chi2 = ((yyteff-teff)/(yyteff*0.01))^2 ;; this assumes a 1% error in the YY tracks

if keyword_set(debug) or keyword_set(psname) then begin
   mydevice=!d.name

   ;; prepare the plotting device
   if keyword_set(psname) then begin
      ;; astrobetter.com tip on making pretty IDL plots
      set_plot, 'PS'
      aspect_ratio=1
      xsize=10.5/1.5
      ysize=xsize/aspect_ratio
      !p.font=0
      device, filename=psname, /color, bits=24,/encapsulated
      device, xsize=xsize,ysize=ysize
      loadct, 39, /silent
      red = 254
      symsize = 0.33
      xtitle=textoidl('T_{eff}')
      ytitle=textoidl('log g_*')
   endif else begin
      red = '0000ff'x
      symsize = 1
      device,window_state=win_state
      if win_state[29] eq 1 then wset, 29 $
      else window, 29, retain=2
      xtitle='T_eff'
      ytitle='log g'
   endelse

   ;; make a publication-ready plot of the YY track -- Teff vs logg
   rstar = yytrack[1,*]
   loggplottrack =  alog10(mstar/(rstar^2)*gravitysun)
   teffplottrack = yytrack[0,*]
   loggplot =  alog10(mstar/(yyrstar^2)*gravitysun)

   xmin=max(teffplottrack,min=xmax) ;; plot range backwards
   ymin = min([loggplot,3,5],max=ymax)

   plot, teffplottrack, loggplottrack,xtitle=xtitle,ytitle=ytitle, xrange=[xmin,xmax], yrange=[ymin,ymax]
   plotsym,0,/fill
   oplot, [teff], [loggplot], psym=8,symsize=0.5 ;; the model point

   if keyword_set(psname) then device, /close
   set_plot, mydevice

endif

return, chi2

stop
end
