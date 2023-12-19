;+
; NAME:
;   STROM_CONV
; PURPOSE:
;    Translates stromgren color combinations from the catalog to
;    individual uvby magnitudes
; Modification 
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
;
;-
function strom_conv,V,sigV,by,sigby,m1,sigm1,c1,sigc1,silent=silent,useticav=useticav

if n_params() lt 8 then begin
  print,'syntax: result=strom_conv(V,sigV,by,sigby,m1,sigm1,c1,sigc1)'
  retall
endif

y = V
b = V + by 
u = V + 3*by + 2*m1 + c1
v = V + 2*by + m1

sigy = sigV
sigb = sqrt(sigV^2 + sigby^2)
sigu = sqrt(sigV^2 + (3*sigby)^2 + (2*sigm1)^2 + sigc1^2)
sigv = sqrt(sigV^2 + (2*sigby)^2 + sigm1^2)

if not keyword_set(silent) then begin
  print,'uvby:',u,v,b,y
  print,'sig_uvby:',sigu,sigv,sigb,sigy
endif

return,[u,sigu,v,sigv,b,sigb,y,sigy]

end

;+
; NAME:
;   MKTICSED
;
; PURPOSE: 
;   Create the SED input file for EXOFASTv2 based on TICv8.2. It also
;   creates a minimal prior file that includes Gaia DR3 parallax,
;   upper limits from Galactic dust maps, and starting values for the
;   star from TICv8.2.
;
; CALLING SEQUENCE:
;   mkticsed, ticid [, PRIORFILE=, SEDFILE=, /FRANCE, RA=, DEC=]
; 
; INPUTS:
;
;  TICID - A string representing the numerical portion of the TICID
;          for your target. This can be retrived from EXOFOP.
;          https://exofop.ipac.caltech.edu/tess/
;
; OPTIONAL INPUTS:
;
;  PRIORFILE - The name of the output prior file. By default, it will
;              be TICID.priors
;  SEDFILE   - The name of the output SED file. By default, it will be
;              TICID.sed 
;  RA/DEC    - The right ascension and declination of the object, in
;              decimal degrees. If RA/DEC are provided, we find the
;              TICID closest to the provided RA/Dec. This can be
;              problematic for high proper motion stars or stars in
;              crowded fields, so it is strongly recommended to use
;              the TICID instead.
;
; OPTIONAL KEYWORDS:
;
;  FRANCE    - Use the FRANCE mirror when querying vizier. By default,
;              we use the CfA mirror.
;  
; EXAMPLES:
;
;  ;; make an sed and prior file for WASP-4:
;  mkticsed, '402026209', priorfile='wasp4.priors', sedfile='wasp4.sed'
;
; MODIFICATION HISTORY
; 
;  2023/02 -- Documented.
;-

pro mkticsed, ticid, priorfile=priorfile, sedfile=sedfile, gaiaspfile=gaiaspfile, france=france, ra=ra, dec=dec, $
              apass=apass, galax=galex, ucac=ucac, merm=merm, stromgren=stromgren, kepler=kepler, tycho=tycho

if !version.os_family eq 'Windows' then $
   message,'This program relies on queryvizier, which is not supported in Windows'

!except=0

;; for use without a license
if lmgr(/vm) or lmgr(/runtime) then begin
   par = command_line_args(count=numargs)
   
   ticid = strtrim(par[0],2)
   for i=0L, numargs-1 do begin
      if strpos(par[i],'=') ne -1 then begin
         entries = strsplit(par[i],'=',/extract)
         if strupcase(entries[0]) eq 'PRIORFILE' then priorfile = strtrim(entries[1],2)
         if strupcase(entries[0]) eq 'SEDFILE' then sedfile = strtrim(entries[1],2)
         if strupcase(entries[0]) eq 'FRANCE' then france = long(entries[1])
         if strupcase(entries[0]) eq 'RA' then ra = double(entries[1])
         if strupcase(entries[0]) eq 'DEC' then dec = double(entries[1])
      endif
   endfor
endif

if keyword_set(france) then cfa = 0B $
else cfa = 1B

;; match RA/Dec to TICv8 (closest)
if n_elements(ra) ne 0 and n_elements(dec) ne 0 then begin
   print, 'WARNING: querying by RA/Dec is less robust than querying by TIC ID and may lead to misidentification'   
;   qtic = exofast_queryvizier('IV/38/tic',[ra,dec],2d0, /silent, /all, cfa=cfa)
   qtic = exofast_queryvizier('IV/39/tic82',[ra,dec],2d0, /silent, /all, cfa=cfa)
   if (size(qtic))[2] ne 8 then begin
      print, 'No match to ' + string(ra,dec,format='(f0.8,",",f0.8)')
      return
   endif
   junk = min(qtic._r,ndx)
   ticid = qtic[ndx].tic
   print, 'Matching TIC ID is ' + strtrim(ticid,2)
endif

if n_elements(ticid) eq 0 then message, 'TICID is required'
if n_elements(priorfile) eq 0 then priorfile = ticid + '.priors'
if n_elements(sedfile) eq 0 then sedfile = ticid + '.sed'
if n_elements(gaiaspfile) eq 0 then gaiaspfile = ticid + '.Gaia.sed'

dist = 120d0
galdist=dist

;; query TICv8 for the TIC ID
if strtrim(long(ticid),2) eq ticid then begin
;   qtic = Exofast_Queryvizier('IV/38/tic','TIC ' + strtrim(ticid,2),/allcolumns,cfa=cfa)
   qtic = Exofast_Queryvizier('IV/39/tic82','TIC ' + strtrim(ticid,2),/allcolumns,cfa=cfa)
endif else begin
   ;; query by supplied name (less robust)
   print, 'WARNING: querying by name is less robust than querying by TIC ID and may lead to misidentification'
;   qtic = Exofast_Queryvizier('IV/38/tic',ticid,2d0,/allcolumns,cfa=cfa)
   qtic = Exofast_Queryvizier('IV/39/tic82',ticid,2d0,/allcolumns,cfa=cfa)
   qgaia = Exofast_Queryvizier('I/345/gaia2',ticid,2d0,/allcolumns,cfa=cfa)   

   sorted = sort(qgaia.gmag)

   nstars = n_elements(qgaia)
   if strpos(ticid,'B') eq (strlen(ticid)-1) and nstars ge 2 then begin
      gaiaid = qgaia[sorted[1]].source
   endif else gaiaid = qgaia[sorted[0]].source
   match = where(qtic.gaia eq gaiaid)
   if match[0] eq -1 then message, 'no matching star found. try using the TIC ID directly'
   qtic = qtic[match]

   ticid = strtrim(qtic.tic,2)
   print, 'Matching TIC ID is ' + strtrim(ticid,2)
endelse

if (size(qtic))[2] ne 8 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

match = (where(qtic.tic eq ticid))[0]
if match eq -1 then begin
   print, 'No match to ' + strtrim(ticid,2)
   return
endif

;; prior file
openw, priorlun, priorfile, /get_lun
printf, priorlun, '#### TICv8.2 ####'

;; warn user if this object is a duplicate
if qtic[match].disp eq 'SPLIT' or qtic[match].disp eq 'DUPLICATE' then begin
   printf, priorlun, '# WARNING: disposition in TICv8.2 is ' + qtic[match].disp
   duplicate_id = qtic[match].m_tic
   ;; do I need to use both IDs? Maybe that's what split means?
   ;; need to read and understand this paper: https://arxiv.org/pdf/2108.04778.pdf
   if duplicate_id ne '-1' then begin
      printf, priorlun, '# WARNING: Using TIC ID of duplicate (' + strtrim(duplicate_id,2) + ')'
      match = (where(qtic.tic eq duplicate_id))[0]
   endif
endif

qtic = qtic[match]
star = [qtic.raj2000,qtic.dej2000]

if finite(qtic.mass) and finite(qtic.rad) and finite(qtic.teff) then begin
   ;; require all three. starting hybrid
   ;; values is less robust than starting at solar
   printf, priorlun, qtic.mass, format='("mstar",x,f0.2)'
   printf, priorlun, qtic.rad, format='("rstar",x,f0.2)'
   printf, priorlun, qtic.teff, format='("teff",x,i5)'
endif

;; use the TICv8 prior on metallicity if available
feh = qtic._m_h_
;; this order preserves NaN, sets floor of 0.08
ufeh = 0.08d0 > qtic.e__m_h_ 
if finite(feh) and finite(ufeh) then begin
   printf, priorlun, feh, ufeh, format='("feh",x,f0.5,x,f0.5)'
endif

;; extinction prior 
;; Use TICv8 Gaussian prior if availble
;; upper limit from the Schlegel dust map if not
av = qtic.e_b_v*3.1d0
uav = (0.02d0 > qtic.s_e_b_v)*3.1d0
if finite(av) and finite(uav) and keyword_set(useticav) then begin
   printf, priorlun, av, uav, format='("av",x,f0.5,x,f0.5)'
   printf, priorlun, '##############'
endif else begin
   printf, priorlun, '##############'
   junk = getavprior(ra=qtic.raj2000, dec=qtic.dej2000, line=line)
   printf, priorlun, line
endelse

;if finite(qtic.teff) then printf, priorlun, qtic.teff, format='("teff",x,i5)'

;; open the SED file for writing
fmt = '(a13,x,f9.6,x,f0.6,x,f0.6)'
openw, lun, sedfile, /get_lun
printf, lun, '# bandname magnitude used_errors catalog_errors'

;; warn user if this object is a duplicate
if qtic.disp eq 'SPLIT' or qtic.disp eq 'DUPLICATE' then begin
   printf, lun, '# WARNING: disposition in TICv8.2 is ' + qtic.disp
endif

;; warn the user if this target appears in the Washington Double Star catalog
qwds = Exofast_Queryvizier('B/wds/wds','TIC ' + ticid,2d0,/allcolumns,cfa=cfa)
if (size(qwds))[2] eq 8 then begin
   printf, lun, '# WARNING: this star appears in the Washington Double Star catalog'
   printf, lun, '# WARNING: this star is a visual binary with separation ' + strtrim(qwds.sep2,2)
   printf, lun, '# WARNING: unresolved photometry will bias the SED and the whole fit'
endif

;; use the Gaia ID to query the Gaia catalog
gaiaid = qtic.gaia
qgaia=Exofast_Queryvizier('I/345/gaia2',star,dist/60.,/silent,cfa=cfa,/all)
dr2str = ''
if (size(qgaia))[2] eq 8 then begin
   match = (where(qgaia.source eq gaiaid))[0]
   if match ne -1 then begin
      qgaia = qgaia[match]
      
      if finite(qgaia.plx) and finite(qgaia.e_plx) and qgaia.plx gt 0d0 then begin
         ;; gaia parallax prior 
         ;; with Stassun & Torres correction
;         printf, priorlun, qgaia.plx + 0.082d0, sqrt(qgaia.e_plx^2 + 0.033d0^2), format='("parallax",x,f0.5,x,f0.5)'
         ;; with Lindegren+ 2018 correction
         k = 1.08d0
         if qgaia.gmag le 13d0 then begin
            sigma_s = 0.021d0
         endif else begin
            sigma_s = 0.043d0
         endelse
         printf, priorlun, "# NOTE: the Gaia DR2 parallax (" + strtrim(qgaia.plx,2) + ") and uncertainty (" + strtrim(qgaia.e_plx,2) + ") has been corrected as prescribed in Lindegren+ (2018)"
         corrected_plx = qgaia.plx + 0.030d0
         if corrected_plx gt 0d0 then begin
            dr2str = string(corrected_plx, sqrt((k*qgaia.e_plx)^2 + sigma_s^2), format='("parallax",x,f0.5,x,f0.5)')
            printf, priorlun, '#' + dr2str
         endif else begin
            printf, priorlun, '# WARNING: Negative parallax is not allowed. This is the corrected (but disallowed) parallax'
            printf, priorlun, corrected_plx, sqrt((k*qgaia.e_plx)^2 + sigma_s^2), format='("#parallax",x,f0.5,x,f0.5)'
            upperlimit = corrected_plx + 3d0*sqrt((k*qgaia.e_plx)^2 + sigma_s^2)
            if upperlimit gt 1d-3 then begin
               printf, priorlun, '# Applying a 3-sigma upper limit instead'
               printf, priorlun, upperlimit, format='("#parallax 0.001 -1 0 ",f0.5)'
            endif else begin
               printf, priorlun, '# 3-sigma upper limit is still negative, ignoring parallax'
            endelse
         endelse
      endif      

      if qgaia.gmag gt -9 and finite(qgaia.e_gmag) and (qgaia.e_gmag lt 1d0) then printf,     lun,'#        Gaia',qgaia.gmag,max([0.02d,qgaia.e_gmag]),qgaia.e_gmag, format=fmt
      if qgaia.bpmag gt -9 and finite(qgaia.e_bpmag) and (qgaia.e_bpmag lt 1d0) then printf,  lun,'#      GaiaBP',qgaia.bpmag,max([0.02d,qgaia.e_bpmag]),qgaia.e_bpmag, format=fmt
      if qgaia.rpmag gt -9 and finite(qgaia.e_rpmag)  and (qgaia.e_rpmag lt 1d0) then printf, lun,'#      GaiaRP',qgaia.rpmag,max([0.02d,qgaia.e_rpmag]),qgaia.e_rpmag, format=fmt      

   endif
endif

;; EDR3 
;qgaia3=Exofast_Queryvizier('I/350/gaiaedr3',star,dist/60.,/silent,cfa=cfa,/all)
;; DR3
qgaia3=Exofast_Queryvizier('I/355/gaiadr3',star,dist/60.,/silent,cfa=cfa,/all)
if (size(qgaia3))[2] eq 8 then begin
   match = (where(qgaia3.source eq gaiaid))[0]
   if match ne -1 then begin
      qgaia3 = qgaia3[match]
      
      if finite(qgaia3.plx) and finite(qgaia3.e_plx) and qgaia3.plx gt 0d0 then begin

         phot_g_mean_mag = qgaia3.gmag 
         nu_eff_used_in_astrometry = qgaia3.nueff
         pseudocolor = qgaia3.pscol
         ecl_lat = qgaia3.elat
         astrometric_params_solved = qgaia3.solved

         ;; is the EDR3 error below the floor? If so, round up
;         if qgaia3.e_plx lt 0.03d0 then begin
;            printf, priorlun, "# NOTE: the Gaia EDR3 parallax error (" + strtrim(qgaia3.e_plx,2) + ") has been rounded up to 30 uas to account for systematic floors described in Lindegren+ 2021"
;            uplx = 0.03d0
;         endif else uplx = qgaai3.e_plx

         uplx = sqrt(qgaia3.e_plx^2 + 0.01d0^2)

         ;; is it within range of the Lindegren+ 2021 prescription?
         if ( (astrometric_params_solved eq 31 and nu_eff_used_in_astrometry ge 1.1d0 and nu_eff_used_in_astrometry le 1.9d0) or $
              (astrometric_params_solved eq 95 and pseudocolor ge 1.24d0 and pseudocolor le 1.72d0)) and $
            phot_g_mean_mag ge 6d0 and phot_g_mean_mag le 21d0 then begin
            zpt = get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolor, ecl_lat, astrometric_params_solved)
            printf, priorlun, "# NOTE: the Gaia DR3 parallax (" + strtrim(qgaia3.plx,2) + ") has been corrected by subtracting " + strtrim(zpt,2) + " mas as prescribed in Lindegren+ 2021"
            printf, priorlun, "# NOTE: the Gaia DR3 parallax uncertainty (" + strtrim(qgaia3.e_plx,2) + ") has been added in quadrature with 0.01 to account for remaining systematic residuals"
            printf, priorlun, qgaia3.plx-zpt, uplx, format='("parallax",x,f0.5,x,f0.5)'            
         endif else begin
            printf, priorlun, "# NOTE: the Gaia DR3 parallax could not be corrected and is raw from the catalog"
            printf, priorlun, qgaia3.plx, uplx, format='("parallax",x,f0.5,x,f0.5)'
         endelse
      endif else if dr2str ne '' then begin
         printf, priorlun, "# DR3 parallax unavailable. Using DR2 parallax"
         printf, priorlun, dr2str
      endif

      if qgaia3.gmag gt -9 and finite(qgaia3.e_gmag) and (qgaia3.e_gmag lt 1d0)  then printf, lun,'Gaia_G_EDR3',qgaia3.gmag,max([0.02d,qgaia3.e_gmag]),qgaia3.e_gmag, format=fmt
      if qgaia3.bpmag gt -9 and finite(qgaia3.e_bpmag) and (qgaia3.e_bpmag lt 1d0) then printf, lun,'Gaia_BP_EDR3',qgaia3.bpmag,max([0.02d,qgaia3.e_bpmag]),qgaia3.e_bpmag, format=fmt
      if qgaia3.rpmag gt -9 and finite(qgaia3.e_rpmag) and (qgaia3.e_rpmag lt 1d0) then printf, lun,'Gaia_RP_EDR3',qgaia3.rpmag,max([0.02d,qgaia3.e_rpmag]),qgaia3.e_rpmag, format=fmt      
   endif else if dr2str ne '' then begin
      printf, priorlun, "# DR3 parallax unavailable. Using DR2 parallax"
      printf, priorlun, dr2str
   endif
endif else if dr2str ne '' then begin
   printf, priorlun, "# DR3 parallax unavailable. Using DR2 parallax"
   printf, priorlun, dr2str
endif

;; Gaia spectrophotometry
qgaiasp=Exofast_Queryvizier('I/355/xpsample',star,dist/60.,/silent,cfa=cfa,/all)
if (size(qgaiasp))[2] eq 8 then begin
   match = where(qgaiasp.source eq gaiaid)
   if match[0] ne -1 then begin
      ;; Gaia lambda in nm, Gaia flux in W/m^2/Hz
      ;; SED plot in microns, erg/s/cm^2
      ;; interpolate onto atmosphere wavelength scale now?
      flux = qgaiasp[match].flux*1d3*qgaiasp[match].lambda ;; W/m^2/Hz -> erg/s/cm^2
      fluxerr = qgaiasp[match].e_flux*1d3*qgaiasp[match].lambda ;; W/m^2/Hz -> erg/s/cm^2
      lambda = qgaiasp[match].lambda/1d3 ;; nm -> um
      exofast_forprint, lambda, flux, fluxerr, $
                        textout=gaiaspfile,/nocomment,format='(f5.3,x,e12.6,x,e12.6)'
   endif
endif

;; use the 2MASS ID to query the 2MASS catalog
tmassid = qtic._2mass
q2mass=Exofast_Queryvizier('II/246/out',star,dist/60.,/silent,cfa=cfa)
if (size(q2mass))[2] eq 8 then begin
   match = (where(q2mass._2mass eq tmassid))[0]
   if match ne -1 then begin
      q2mass = q2mass[match]
      if q2mass.Jmag gt -9 and (q2mass.e_Jmag lt 1d0) then printf, lun,'J2M',q2mass.Jmag,max([0.02d,q2mass.e_Jmag]),q2mass.e_Jmag, format=fmt
      if q2mass.Hmag gt -9 and (q2mass.e_Hmag lt 1d0) then printf, lun,'H2M',q2mass.Hmag,max([0.02d,q2mass.e_Hmag]),q2mass.e_Hmag, format=fmt
      if q2mass.Kmag gt -9 and (q2mass.e_Kmag lt 1d0) then begin
         printf, lun,'K2M',q2mass.Kmag,max([0.02d,q2mass.e_Kmag]),q2mass.e_Kmag, format=fmt
         ;; include the 
         printf, priorlun,'# Apparent 2MASS K magnitude for the Mann relation'
         printf, priorlun,'appks',q2mass.Kmag,max([0.02d,q2mass.e_Kmag]), format='(a5,x,f9.6,x,f0.6)'
      endif
   endif
endif

;; use the WISE ID to query the wise catalog
wiseid = qtic.wisea
qwise=Exofast_Queryvizier('II/328/allwise',star,dist/60.,/silent,cfa=cfa)
if (size(qwise))[2] eq 8 then begin
   match = (where(qwise.allwise eq wiseid))[0]
   if match ne -1 then begin
      qwise = qwise[match]
      if qwise.w1mag gt -9 and finite(qwise.e_w1mag) and (qwise.e_w1mag lt 1d0) then printf, lun,'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag]),qwise.e_w1mag, format=fmt
      if qwise.w2mag gt -9 and finite(qwise.e_w2mag) and (qwise.e_w2mag lt 1d0) then printf, lun,'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag]),qwise.e_w2mag, format=fmt
      if qwise.w3mag gt -9 and finite(qwise.e_w3mag) and (qwise.e_w3mag lt 1d0) then printf, lun,'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag]),qwise.e_w3mag, format=fmt
      if qwise.w4mag gt -9 and finite(qwise.e_w4mag) and (qwise.e_w4mag lt 1d0) then printf, lun,'WISE4',qwise.w4mag,max([0.10d,qwise.e_w4mag]),qwise.e_w4mag, format=fmt
   endif else begin
      print, 'No match in WISE by ID in TICv8.2'
      mindmag = min(abs(q2mass.Kmag-qwise.w1mag),match)
      sep = angsep(qtic.raj2000*!dpi/180d0,qtic.dej2000*!dpi/180d0,qwise.raj2000*!dpi/180d0, qwise.dej2000*!dpi/180d0)*180d0/!dpi*3600d0 ;; arcsec
      if mindmag lt 0.5 and sep[match] lt 30d0 then begin
         qwise = qwise[match]
         printf, lun, '# No match in WISE by ID in TICv8.2, matched by K-WISE1 mag (' + strtrim(mindmag,2) + ') and separation (' + strtrim(sep[match],2) + '")'
         if qwise.w1mag gt -9 and finite(qwise.e_w1mag) and (qwise.e_w1mag lt 1d0) then printf, lun,'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag]),qwise.e_w1mag, format=fmt
         if qwise.w2mag gt -9 and finite(qwise.e_w2mag) and (qwise.e_w2mag lt 1d0) then printf, lun,'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag]),qwise.e_w2mag, format=fmt
         if qwise.w3mag gt -9 and finite(qwise.e_w3mag) and (qwise.e_w3mag lt 1d0) then printf, lun,'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag]),qwise.e_w3mag, format=fmt
         if qwise.w4mag gt -9 and finite(qwise.e_w4mag) and (qwise.e_w4mag lt 1d0) then printf, lun,'WISE4',qwise.w4mag,max([0.10d,qwise.e_w4mag]),qwise.e_w4mag, format=fmt
      endif
   endelse
endif

;; use the Tycho ID to query the Stromgren catalog to get a metalicity prior
if ~finite(feh) or ~finite(ufeh) then begin
   tycid = qtic.tyc
   qpaunzen15=Exofast_Queryvizier('J/A+A/580/A23/catalog',star,dist/60.,/silent,/all,cfa=cfa)
   if (size(qpaunzen15))[2] eq 8 then begin
      paunzentycid = string(qpaunzen15.tyc1,qpaunzen15.tyc2,qpaunzen15.tyc3,format='(i04,"-",i05,"-",i01)')
      match = (where(paunzentycid eq tycid))[0]
      if match ne -1 then begin
         qpaunzen15 = qpaunzen15[match]
         
         ;; don't use these for the SED (?)
         if 0 then begin
            ubvymags=strom_conv(qpaunzen15.vmag,max([0.01d,qpaunzen15.e_vmag]),$
                                qpaunzen15.b_y,max([0.02d,qpaunzen15.e_b_y]),$
                                qpaunzen15.m1,max([0.02d,qpaunzen15.e_m1]),$
                                qpaunzen15.c1,max([0.02d,qpaunzen15.e_c1]),/silent)
            printf, lun, '# Stromgren photometry, Paunzen, 2015'
            printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
            
            if (ubvymags(0) gt -9) and (ubvymags(1) lt 1d0) then printf, lun,'uStr',ubvymags(0),max([0.02d,ubvymags(1)]), format=fmt
            if (ubvymags(2) gt -9) and (ubvymags(3) lt 1d0) then printf, lun,'vStr',ubvymags(2),max([0.02d,ubvymags(3)]), format=fmt
            if (ubvymags(4) gt -9) and (ubvymags(5) lt 1d0) then printf, lun,'bStr',ubvymags(4),max([0.02d,ubvymags(5)]), format=fmt
            if (ubvymags(6) gt -9) and (ubvymags(7) lt 1d0) then printf, lun,'yStr',ubvymags(6),max([0.02d,ubvymags(7)]), format=fmt
         endif
         
         b_y = qpaunzen15.b_y
         m1 = qpaunzen15.m1
         c1 = qpaunzen15.c1
         
         ;; metalicity prior from Cassegrande+ 2011 (solar neighborhood)
         if b_y gt 0.23d0 and b_y lt 0.63d0 and $
            m1 gt 0.05d0 and m1 le 0.68d0 and $
            c1 gt 0.13d0 and c1 le 0.60d0 then begin
            
            ;; Cassegrande+ 2011, eq 2
            feh = 3.927d0*alog10(m1) - 14.459d0*m1^3 - 5.394d0*b_y*alog10(m1) + $
                  36.069d0*b_y*m1^3 + 3.537d0*c1*alog10(m1) - $
                  3.500d0*m1^3*c1 + 11.034d0*b_y - 22.780d0*b_y^2 + $
                  10.684d0*c1 - 6.759d0*c1^2 - 1.548d0
            ufeh = 0.10d0
            
            printf, lun, '# Stromgren photometry, Paunzen, 2015'
            printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
            printf, priorlun, '# Cassegrande+ 2011, eq 2'
            printf, priorlun, feh, ufeh, format='("feh",x,f0.5,x,f0.5)'
         endif else if b_y gt 0.43d0 and b_y lt 0.63d0 and $
            m1 gt 0.07d0 and m1 le 0.68d0 and $
            c1 gt 0.16d0 and c1 le 0.49d0 then begin
            
            ;; Cassegrande+ 2011, eq 3
            feh = -0.116d0*c1 - 1.624d0*c1^2 + 8.955d0*c1*b_y + $
                  42.008d0*b_y - 99.596d0*b_y^2 + 64.245d0*b_y^3 + $
                  8.928d0*c1*m1 + 17.275d0*m1 - 48.106d0*m1^2 + $
                  45.802d0*m1^3 - 8.467d0
            ufeh = 0.12d0
            printf, lun, '# Stromgren photometry, Paunzen, 2015'
            printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P' 
            printf, priorlun, '# Cassegrande+ 2011, eq 3'
            printf, priorlun, feh, ufeh, format='("feh",x,f0.5,x,f0.5)'
         endif      
         
      endif
   endif
endif

;; we can do better than -5 < [Fe/H]_0 < +0.5
if ~finite(feh) or ~finite(ufeh) then begin
;   printf, priorlun, '# Cassegrande+ 2011, Table 1'
;   feh = -0.06d0
;   ufeh = 0.25d0
   
   printf, priorlun, '# wide Gaussian prior'
   feh = 0d0
   ufeh = 1d0
   printf, priorlun, feh, ufeh, format='("feh",x,f0.5,x,f0.5)'
endif

;; get GALEX -- by default, skip; models are unreliable here
qgalex=Exofast_Queryvizier('II/312/ais',star,galdist/60d0,/silent,/all,cfa=cfa)
if long(tag_exist(qgalex,'fuv',/quiet)) ne 0L then begin
   printf, lun, '# Galex DR5, Bianchi+ 2011'
   printf, lun, '# http://adsabs.harvard.edu/abs/2011Ap%26SS.335..161B'
   printf, lun, '# Atmospheric models are generally untrustworthy here, and including UV photometry may bias the fit'
   printf, lun, '# Matching is done by nearest neighbor with ~10% failure rate'
   if n_elements(qgalex) gt 1 then begin
      print,'Warning: More than 1 GALEX source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 GALEX source found; using nearest one only.'
      junk = min(qgalex._r,m)
      qgalex=qgalex[m]
   endif
   if qgalex.fuv gt -99 and finite(qgalex.e_fuv) then begin
      line = string('galFUV',qgalex.fuv,max([0.1d,qgalex.e_fuv]),qgalex.e_fuv, format=fmt) 
      if keyword_set(galex) then printf, lun, line $
      else printf, lun, '#' + line
   endif
   if qgalex.nuv gt -99 and finite(qgalex.e_nuv) then begin
      line = string('galNUV',qgalex.nuv,max([0.1d,qgalex.e_nuv]),qgalex.e_nuv, format=fmt)
      if keyword_set(galex) then printf, lun,line $
      else printf, lun, '#' + line
   endif
endif else qgalex={fuv_6:-99.,nuv_6:-99.}

; Tycho-2
qtyc2=Exofast_Queryvizier('I/259/TYC2',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qtyc2,'BTMAG',/quiet)) ne 0L then begin
   printf, lun, '# Tycho catalog, Hoeg+ 2000'
   printf, lun, '# http://adsabs.harvard.edu/abs/2000A%26A...355L..27H'
   printf, lun, '# Matching is done by nearest neighbor with ~10% failure rate'
   if n_elements(qtyc2) gt 1 then begin
      print,'Warning: More than 1 Tycho-2 source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Tycho-2 source found; using nearest one only.'
      junk = min(qtyc2._r,m)
      qtyc2=qtyc2[m]
   endif
   if keyword_set(tycho) then comment = '' else comment = '#          '
   if qtyc2.btmag gt -9 and finite(qtyc2.e_btmag) then printf, lun,comment+'BT',qtyc2.btmag,max([0.02d,qtyc2.e_btmag]),qtyc2.e_btmag, format=fmt
   if qtyc2.vtmag gt -9 and finite(qtyc2.e_vtmag) then printf, lun,comment+'VT',qtyc2.vtmag,max([0.02d,qtyc2.e_vtmag]),qtyc2.e_vtmag, format=fmt
endif else qtyc2={btmag:-99.,vtmag:-99.}

; UCAC4 
qucac4=Exofast_Queryvizier('UCAC4',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qucac4,'bmag',/quiet)) ne 0L then begin
   printf, lun,'# UCAC4, Zacharias+, 2012'
   printf, lun,'# http://adsabs.harvard.edu/abs/2012yCat.1322....0Z'
   printf, lun,'# Matching is done by nearest neighbor with ~10% failure rate'
   if n_elements(qucac4) gt 1 then begin
      print,'Warning: More than 1 UCAC-4 source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 UCAC-4 source found; using nearest one only.'
      junk = min(qucac4._r,m)
      qucac4=qucac4[m]
   endif
   if keyword_set(apass) then comment1 = '' else comment1 = '#           ' 
   if keyword_set(apass) then comment2 = '' else comment2 = '#       ' 

   printf, lun, '# APASS DR6 (via UCAC4), Henden+ 2016'
   printf, lun, '# http://adsabs.harvard.edu/abs/2016yCat.2336....0H'
   if qucac4.bmag ne qtyc2.btmag and qucac4.bmag gt -9 and qucac4.e_bmag ne 99 then printf, lun,comment1+'B',qucac4.bmag,max([0.02d,qucac4.e_bmag*0.01d]),qucac4.e_bmag*0.01d, format=fmt
   if qucac4.vmag ne qtyc2.vtmag and qucac4.bmag gt -9 and qucac4.e_vmag ne 99 then printf, lun,comment1+'V',qucac4.vmag,max([0.02d,qucac4.e_vmag*0.01d]),qucac4.e_vmag*0.01d, format=fmt
   if qucac4.gmag gt -9 then printf, lun,comment2+'gSDSS',qucac4.gmag,max([0.02d,qucac4.e_gmag*0.01d]),qucac4.e_gmag*0.01d, format=fmt
   if qucac4.rmag gt -9 then printf, lun,comment2+'rSDSS',qucac4.rmag,max([0.02d,qucac4.e_rmag*0.01d]),qucac4.e_rmag*0.01d, format=fmt
   if qucac4.imag gt -9 then printf, lun,comment2+'iSDSS',qucac4.imag,max([0.02d,qucac4.e_imag*0.01d]),qucac4.e_imag*0.01d, format=fmt

;   if qucac4.Jmag gt -9  or qucac4.Hmag gt -9 or qucac4.Kmag gt -9 then begin
;      printf, lun, '# 2MASS (via UCAC4), Cutri+ 2003'
;      printf, lun, '# http://adsabs.harvard.edu/abs/2003yCat.2246....0C'   
;   endif
;   if qucac4.Jmag gt -9 then printf, lun,'J2M',qucac4.Jmag,max([0.02d,qucac4.e_Jmag]),qucac4.e_Jmag, format=fmt
;   if qucac4.Hmag gt -9 then printf, lun,'H2M',qucac4.Hmag,max([0.02d,qucac4.e_Hmag]),qucac4.e_Hmag, format=fmt
;   if qucac4.Kmag gt -9 then printf, lun,'K2M',qucac4.Kmag,max([0.02d,qucac4.e_Kmag]),qucac4.e_Kmag, format=fmt
endif; else qucac4={Jmag:-99.,Hmag:-99.,Kmag:-99.}


; Paunzen, 2015 Stromgren
qpaunzen15=Exofast_Queryvizier('J/A+A/580/A23/catalog',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qpaunzen15,'vmag',/quiet)) ne 0L then begin
   if n_elements(qpaunzen15) gt 1 then begin
      print,'Warning: More than 1 Paunzen source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Paunzen source found; using nearest one only.'
      junk = min(qpaunzen15._r,m)
      qpaunzen15=qpaunzen15[m]
   endif
   ubvymags=strom_conv(qpaunzen15.vmag,max([0.01d,qpaunzen15.e_vmag]),$
                       qpaunzen15.b_y,max([0.02d,qpaunzen15.e_b_y]),$
                       qpaunzen15.m1,max([0.02d,qpaunzen15.e_m1]),$
                       qpaunzen15.c1,max([0.02d,qpaunzen15.e_c1]),/silent)
   printf, lun, '# Stromgren photometry, Paunzen, 2015'
   printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
   printf, lun, '# Matching is done by nearest neighbor with ~10% failure rate'
   if keyword_set(stromgren) then comment = '' else comment = '#        '

   if ubvymags(0) gt -9 then printf, lun,comment+'uStr',ubvymags(0),max([0.02d,ubvymags(1)]), format=fmt
   if ubvymags(2) gt -9 then printf, lun,comment+'vStr',ubvymags(2),max([0.02d,ubvymags(3)]), format=fmt
   if ubvymags(4) gt -9 then printf, lun,comment+'bStr',ubvymags(4),max([0.02d,ubvymags(5)]), format=fmt
   if ubvymags(6) gt -9 then printf, lun,comment+'yStr',ubvymags(6),max([0.02d,ubvymags(7)]), format=fmt
endif else qpaunzen15={vmag:-99.}

; get Mermilliod 1991 UBV ;; by default skip (Gaia better)
qmermilliod91=Exofast_Queryvizier('II/168/ubvmeans',star,dist/60.,/silent,/all,cfa=cfa) 
if long(tag_exist(qmermilliod91,'vmag',/quiet)) ne 0L then begin
   if n_elements(qmermilliod91) gt 1 then begin
      print,'Warning: More than 1 Mermilliod source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Mermilliod source found; using nearest one only.'
      junk = min(qmermilliod91._r,m)
      qmermilliod91=qmermilliod91[m]
   endif

   if keyword_set(merm) then comment='' else comment='#           '
   printf, lun, '# Mermilliod, 1994'
   printf, lun, '# http://adsabs.harvard.edu/abs/1994yCat.2193....0M'
   printf, lun, '# Matching is done by nearest neighbor with ~10% failure rate'

   if qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag gt -9 and $
      finite(qmermilliod91.e_u_b) and finite(qmermilliod91.e_b_v) and $
      finite(qmermilliod91.e_vmag) then begin
      printf, lun, comment+'U',qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag,max([0.02d,sqrt(qmermilliod91.e_u_b^2+qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)]), format=fmt
   endif
   
   if qmermilliod91.b_v+qmermilliod91.vmag gt -9 and $
      finite(qmermilliod91.e_b_v) and $
      finite(qmermilliod91.e_vmag) then begin
      printf, lun, comment+'B',qmermilliod91.b_v+qmermilliod91.vmag,max([0.02d,sqrt(qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)]), format=fmt
   endif
   if qmermilliod91.vmag gt -9 and finite(qmermilliod91.e_vmag) then begin
      printf, lun, comment+'V',qmermilliod91.vmag,max([0.02d,qmermilliod91.e_vmag]), format=fmt
   endif
endif else qmermilliod91={vmag:-99.}

; KIS DR2
qkis=Exofast_Queryvizier('J/AJ/144/24/kisdr2',star,dist/60.,/silent,cfa=cfa,/all)
if long(tag_exist(qkis,'KIS',/quiet)) ne 0L then begin
   printf, lun, '# KIS DR2, Greiss+ 2012'
   printf, lun, '# http://adsabs.harvard.edu/abs/2012AJ....144...24G'
   printf, lun, '# Matching is done by nearest neighbor with ~10% failure rate'
   if n_elements(qkis) gt 1 then begin
      print,'Warning: More than 1 KIS source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 KIS source found; using nearest one only.'
      match = crossref(ra, dec, qkis.raj2000, qkis.dej2000)
      qkis = qkis[match]
   endif
   
   if keyword_set(kepler) then comment = '' else comment = '#        '
   
   if qkis.umag gt -9 and finite(qkis.e_umag) then printf, lun,comment+'UKIS',qkis.umag,max([0.02d,qkis.e_umag]),qkis.e_umag, format=fmt
   if qkis.gmag gt -9 and finite(qkis.e_gmag) then printf, lun,comment+'gKIS',qkis.gmag,max([0.02d,qkis.e_gmag]),qkis.e_gmag, format=fmt
   if qkis.rmag gt -9 and finite(qkis.e_rmag) then printf, lun,comment+'rKIS',qkis.rmag,max([0.02d,qkis.e_rmag]),qkis.e_rmag, format=fmt
   if qkis.imag gt -9 and finite(qkis.e_imag) then printf, lun,comment+'iKIS',qkis.imag,max([0.02d,qkis.e_imag]),qkis.e_imag, format=fmt
endif else qkis={umag:-99.,gmag:-99.,rmag:-99.,imag:-99.}



free_lun, priorlun
free_lun, lun

print, 'Successfully retrieved SED and prior information. See ' + priorfile + ' and ' + sedfile

end
