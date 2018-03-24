;; t - teff
;; logg - the logg of the star
;; z  - metalicity
;; wav - wavelength, in ??
;; lamflam -
;; blamflam - 
;; dusty

;; marvels_interp_model, teff, logg, feh

pro exofast_interp_model3d,teff,logg,feh,wav,lamflam,alpha=alpha,dusty=dusty,nextgen=nextgen,cond=cond,kurucz=kurucz,verbose=verbose, interpfiles=interpfiles,logname=logname

errmsg = string(logg, teff, feh, format='("This star (logg=",f0.2,",teff=",f0.2,",[Fe/H]=",f0.2,") is not included in the atmospheric models, likely because it is unphysical. If this is a physical star, turn off SED fitting or it will bias your results")')
if teff lt 800d0 or teff gt 70000d0 then begin
   if keyword_set(verbose) then printandlog, errmsg, logname
   lamflam = !values.d_nan
   return
endif

;; loading files is **SLOW**
;; speed it up by keeping files loaded in this common block
common interp_model, model, specfiles, w1

if n_elements(alpha) eq 0 then alpha = 0.0

nextgen=1
if keyword_set(nextgen) then nextgenpath=filepath('',root_dir=getenv("EXOFAST_PATH"),subdir=['sed','nextgenfin'])

if n_elements(model) eq 0 then begin
   specfiles = file_search(nextgenpath + '*.spec.idl',count=nfiles)
   model = ptrarr(nfiles,/allocate_heap)
   w1 = findgen(24000)/1000+0.1 ;; wavelength scale on which to interpolate
endif

;; select the nearest models
allowedz = [-4d0,-3.5d0,-3d0,-2.5d0,-2d0,-1.5d0,-1d0,-0.5d0,0d0,0.3d0,0.5d0]
if feh lt allowedz[0] or feh gt allowedz[-1] then begin
   if keyword_set(verbose) then printandlog, errmsg, logname
   lamflam = !values.d_nan
   return
endif
zstr = ['-4.0','-3.5','-3.0','-2.5','-2.0','-1.5','-1.0','-0.5','+0.0','+0.3','+0.5']
junk = min(abs(allowedz-feh),zndx)
if allowedz[zndx] gt feh then zndx -= 1

allowedalpha = [0d0,0.2d0,-0.2d0,0.4d0,0.6d0]
alphastr = ['+0.0','+0.2','-0.2','+0.4','+0.6']
nalpha = n_elements(allowedalpha)

allowedlogg = [-0.5d, 0d0, 0.5d0, 1d0, 1.5d0,  2d0, 2.5d0, 3d0, 3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0]
if logg lt allowedlogg[0] or logg gt allowedlogg[-1] then begin
   if keyword_set(verbose) then printandlog, errmsg, logname
   lamflam = !values.d_nan
   return
endif
loggstr = ['-0.5','+0.0','+0.5','+1.0','+1.5','+2.0','+2.5','+3.0','+3.5','+4.0','+4.5','+5.0','+5.5','+6.0']
junk = min(abs(allowedlogg-logg),loggndx)
if allowedlogg[loggndx] gt logg then loggndx -= 1

allowedteff = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, $
               23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, $
               38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, $
               53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, $
               68, 69, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, $
               96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, $
               125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185,$
               190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,$
               310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430,$
               440, 450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560,$
               570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690,$
               700]
teffstr = strtrim(allowedteff,2)
allowedteff *= 100d0
nteff = n_elements(allowedteff)
junk = min(abs(allowedteff-teff),teffndx)
if allowedteff[teffndx] gt teff then teffndx -= 1

if keyword_set(nextgen) or keyword_set(dusty) or keyword_set(cond) then begin

   lamflams = dblarr(24000,2,2,2)
   for i=0, 1 do begin
      for j=0, 1 do begin
         for k=0, 1 do begin
            for l=0, nalpha-1 do begin
               filename = nextgenpath + 'lte' + teffstr[teffndx + i] + loggstr[loggndx+j] + zstr[zndx+k] + alphastr[l]+'.NextGen.spec.idl'
               match = (where(specfiles eq filename))[0]
               if match ne -1 then begin
                  if n_elements(*model[match]) eq 0 then begin
                     restore, specfiles[match]
                     lamflams[*,i,j,k] = lamflam1
                     *model[match] = lamflam1 ;; populate the model variable for next time
                  endif else lamflams[*,i,j,k] = *model[match]
                  break

               endif 
            endfor
            ;; the model file doesn't exist. Unphysical?
            if match[0] eq -1 then begin
               if keyword_set(verbose) then printandlog, errmsg, logname
               lamflam = !values.d_nan
               return
            endif
         endfor
      endfor
   endfor
endif

x_teff = (teff - allowedteff[teffndx])/(allowedteff[teffndx+1] - allowedteff[teffndx])
y_logg = (logg - allowedlogg[loggndx])/(allowedlogg[loggndx+1] - allowedlogg[loggndx])
z_feh  = ( feh -       allowedz[zndx])/(      allowedz[zndx+1] - allowedz[zndx])

lamflam = interpolate(lamflams,x_teff,y_logg,z_feh)

end

