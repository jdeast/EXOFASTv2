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

pro mkticsed, ticid, priorfile=priorfile, sedfile=sedfile

if n_elements(ticid) eq 0 then ticid = '455139555'
if n_elements(priorfile) eq 0 then priorfile = ticid + '.priors'
if n_elements(sedfile) eq 0 then sedfile = ticid + '.sed'

dist = 120d0

;; query TICv8 for the TIC ID
if strtrim(long(ticid),2) eq ticid then begin
   qtic = Exofast_Queryvizier('IV/38/tic','TIC ' + strtrim(ticid,2),/allcolumns,/cfa)
endif else begin
   qtic = Exofast_Queryvizier('IV/38/tic',ticid,2d0,/allcolumns,/cfa)

   qgaia = Exofast_Queryvizier('I/345/gaia2',ticid,2d0,/allcolumns,/cfa)   
   sorted = sort(qgaia.gmag)
   nstars = n_elements(qgaia)
   if strpos(ticid,'B') eq (strlen(ticid)-1) and nstars ge 2 then begin
      gaiaid = qgaia[sorted[1]].source
   endif else gaiaid = qgaia[sorted[0]].source
   match = where(qtic.gaia eq gaiaid)
   if match[0] eq -1 then message, 'no matching star found. try using the TIC ID directly'
   qtic = qtic[match]


; gt -9 and finite(qgaia.e_gmag) then printf, lun,'Gaia',qgaia.gmag,max([0.02d,qgaia.e_gmag]),qgaia.e_gmag, format=fmt
;      if qgaia.bpmag gt -9 and finite(qgaia.e_bpmag) then printf, lun,'GaiaBP',qgaia.bpmag,max([0.02d,qgaia.e_bpmag]),qgaia.e_bpmag, format=fmt
;      if qgaia.rpmag gt -9 and finite(qgaia.e_rpmag) then printf, lun,'GaiaRP',qgaia.rpmag,max([0.02d,qgaia.e_rpmag]),qgaia.e_rpmag, format=fmt      
;      endif

;stop


;gaiaid = qtic.gaia
;qgaia=Exofast_Queryvizier('I/345/gaia2',star,dist/60.,/silent,cfa=cfa,/all)
;if (size(qgaia))[2] eq 8 then begin
;   match = (where(qgaia.source eq gaiaid))[0]
;   if match ne -1 then begin
;      qgaia = qgaia[match]




;   sorted = sort(qtic.tmag)
;   nstars = n_elements(qtic)
;
;   forprint, qtic[sorted].tic, qtic[sorted].tmag
;
;   if strpos(ticid,'B') eq (strlen(ticid)-1) and nstars ge 2 then begin
;      ;; B component -- grab the second brightest thing 
;      qtic = qtic[sorted[1]]
;   endif else qtic = qtic[sorted[0]]
   ticid = strtrim(qtic.tic,2)
print, ticid
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

qtic = qtic[match]
star = [qtic.raj2000,qtic.dej2000]

;; prior file
openw, priorlun, priorfile, /get_lun
printf, priorlun, '#### TICv8 ####'
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
fmt = '(a10,x,f9.6,x,f0.6,x,f0.6)'
openw, lun, sedfile, /get_lun
print, lun, '# bandname magnitude used_errors catalog_errors'

;; use the Gaia ID to query the Gaia catalog
gaiaid = qtic.gaia
qgaia=Exofast_Queryvizier('I/345/gaia2',star,dist/60.,/silent,cfa=cfa,/all)
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
         printf, priorlun, "# NOTE: the Gaia DR2 parallax and uncertainty has been corrected as prescribed in Lindegren+ (2018)"
         printf, priorlun, qgaia.plx + 0.030d0, sqrt((k*qgaia.e_plx)^2 + sigma_s^2), format='("parallax",x,f0.5,x,f0.5)'
      endif      

      if qgaia.gmag gt -9 and finite(qgaia.e_gmag) then printf, lun,'Gaia',qgaia.gmag,max([0.02d,qgaia.e_gmag]),qgaia.e_gmag, format=fmt
      if qgaia.bpmag gt -9 and finite(qgaia.e_bpmag) then printf, lun,'GaiaBP',qgaia.bpmag,max([0.02d,qgaia.e_bpmag]),qgaia.e_bpmag, format=fmt
      if qgaia.rpmag gt -9 and finite(qgaia.e_rpmag) then printf, lun,'GaiaRP',qgaia.rpmag,max([0.02d,qgaia.e_rpmag]),qgaia.e_rpmag, format=fmt      
   endif
endif

;; use the 2MASS ID to query the 2MASS catalog
tmassid = qtic._2mass
q2mass=Exofast_Queryvizier('II/246/out',star,dist/60.,/silent,cfa=cfa)
if (size(q2mass))[2] eq 8 then begin
   match = (where(q2mass._2mass eq tmassid))[0]
   if match ne -1 then begin
      q2mass = q2mass[match]
      if q2mass.Jmag gt -9 then printf, lun,'J2M',q2mass.Jmag,max([0.02d,q2mass.e_Jmag]),q2mass.e_Jmag, format=fmt
      if q2mass.Hmag gt -9 then printf, lun,'H2M',q2mass.Hmag,max([0.02d,q2mass.e_Hmag]),q2mass.e_Hmag, format=fmt
      if q2mass.Kmag gt -9 then printf, lun,'K2M',q2mass.Kmag,max([0.02d,q2mass.e_Kmag]),q2mass.e_Kmag, format=fmt
   endif
endif

;; use the WISE ID to query the wise catalog
wiseid = qtic.wisea
qwise=Exofast_Queryvizier('II/328/allwise',star,dist/60.,/silent,cfa=cfa)
if (size(qwise))[2] eq 8 then begin
   match = (where(qwise.allwise eq wiseid))[0]
   if match ne -1 then begin
      qwise = qwise[match]
      if qwise.w1mag gt -9 and finite(qwise.e_w1mag) then printf, lun,'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag]),qwise.e_w1mag, format=fmt
      if qwise.w2mag gt -9 and finite(qwise.e_w2mag) then printf, lun,'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag]),qwise.e_w2mag, format=fmt
      if qwise.w3mag gt -9 and finite(qwise.e_w3mag) then printf, lun,'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag]),qwise.e_w3mag, format=fmt
      if qwise.w4mag gt -9 and finite(qwise.e_w4mag) then printf, lun,'WISE4',qwise.w4mag,max([0.1d,qwise.e_w4mag]),qwise.e_w4mag, format=fmt
   endif
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
            
            if ubvymags(0) gt -9 then printf, lun,'uStr',ubvymags(0),max([0.02d,ubvymags(1)]), format=fmt
            if ubvymags(2) gt -9 then printf, lun,'vStr',ubvymags(2),max([0.02d,ubvymags(3)]), format=fmt
            if ubvymags(4) gt -9 then printf, lun,'bStr',ubvymags(4),max([0.02d,ubvymags(5)]), format=fmt
            if ubvymags(6) gt -9 then printf, lun,'yStr',ubvymags(6),max([0.02d,ubvymags(7)]), format=fmt
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

free_lun, priorlun
free_lun, lun

end
