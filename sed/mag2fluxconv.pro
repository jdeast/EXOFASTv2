pro	mag2fluxconv,fname,lameff,weff,flux,fluxerr,teff=teff,wiseab=wiseab

; fname - filename with band, mag, magerr
; Returns fluxes and errors in lamFlam units: erg/cm^2/s
; Returns lameff,weff in um
;
; REFS
; Spitzer cookbook unit conversions: 
;   http://www.ipac.caltech.edu/2mass/releases/allsky/faq.html#jansky
; SDSS flux calibration:
;   http://www.sdss.org/dr7/algorithms/fluxcal.html
; Conversion between AB and Vega mags
;   http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
; Tycho flux calibration:
;   http://nsted.ipac.caltech.edu/NStED/docs/parhelp/Photometry.html
; GALEX flux calibration:
;   http://galexgi.gsfc.nasa.gov/docs/galex/Documents/ERO_data_description_2.htm#_Toc58822546
; GALEX fluxes from:
;   http://galex.stsci.edu/GR6/?page=mastform
;   Generally use aper_mag_7, 
;   add 20.08 to NUV instrumental mag,
;   add 18.82 to FUV instrumental mag, 
;   (see http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html)
;   (see also email from Susan Neff)
;   then subtract aperture correction of 0.07 mag, see
;   /home/stassuk/idl/MARVELS/Morrissey_et_al_2007_GALEX_apercorr_fig4.gif
; WISE flux calibration: 
;   http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html#WISEZMA
;   All reported magnitudes are in the Vega System. Need to convert to AB system:
;   add 2.683 to W1
;   add 3.319 to W2
;   add 5.242 to W3
;   add 6.604 to W4
; WISE fluxes from: 
;   http://wise2.ipac.caltech.edu
; 

if n_params() lt 1 then begin
  print,'syntax: mag2fluxconv,fname,lameff,weff,flux,fluxerr,teff=teff'
  retall
endif

if not keyword_set(teff) then teff=5040.
theta=5040./teff

readcol,fname,band,mag,merr,format='a,f,f',/silent,comment='#'
band=strtrim(band,2)
lameff=fltarr(n_elements(band))
weff=fltarr(n_elements(band))
flux=fltarr(n_elements(band))

for i=0,n_elements(band)-1 do begin
  case strn(band(i)) of
     'U':  begin
             lameff(i)=poly(theta,[3476.,162.,86.,-63.]) / 1.e4	; um (ADPS2)
             weff(i)=poly(theta,[612.,346.,-741.,269.]) / 1.e4   ; um (ADPS2)
             lamflamzp = 1.51e-5	; USNO Pub XXV-P1 (1984)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'B':  begin
             lameff(i)=poly(theta,[4336.,201.,235.,-115.]) / 1.e4
             weff(i)=poly(theta,[863.,494.,-833.,192.]) / 1.e4
             lamflamzp = 2.90e-5
             flux(i)=lamflamzp*10^(-0.4*mag(i))
         end
     'V':  begin
             lameff(i)=poly(theta,[5442.,130.,159.,-53.]) / 1.e4
             weff(i)=poly(theta,[827.,167.,-156.,-63.]) / 1.e4
             lamflamzp = 1.98e-5
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'R':  begin
             lameff(i)=poly(theta,[6622.,346.,427.,-202.]) / 1.e4
             weff(i)=poly(theta,[1788.,772.,-361.,-369.]) / 1.e4
             lamflamzp = 1.13e-5
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'RC': begin
             lameff(i)=poly(theta,[6326.,161.,219.,-94.]) / 1.e4
             weff(i)=poly(theta,[1282.,302.,-167.,-147.]) / 1.e4
             lamflamzp = 1.44e-5
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'IC': begin
             lameff(i)=poly(theta,[7814.,64.,79.,-30.]) / 1.e4
             weff(i)=poly(theta,[984.,59.,-3.,-44.]) / 1.e4
             lamflamzp = 9.68e-6
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'J2M': begin
             lameff(i)=poly(theta,[1.23,0.08,0.0,-0.01])
             weff(i)=poly(theta,[0.16,0.08,-0.23,0.11])
             lamflamzp = 1594. * 3e-9 / 1.235^2 * 1.235	; Spitzer cookbook website (Fnu to Flam)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'H2M': begin
             lameff(i)=poly(theta,[1.64,0.05,0.02,-0.02])
             weff(i)=poly(theta,[0.24,0.09,-0.21,0.07])
             lamflamzp = 1024. * 3e-9 / 1.662^2 * 1.662	; Spitzer cookbook website (Fnu to Flam)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'K2M': begin
             lameff(i)=poly(theta,[2.15,0.03,0.02,-0.01])
             weff(i)=poly(theta,[0.25,0.04,-0.06,0.01])
             lamflamzp = 666.7 * 3e-9 / 2.159^2 * 2.159	; Spitzer cookbook website (Fnu to Flam)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'BT': begin
             lameff(i)=4206.4 / 1e4
             weff(i)= 709.8 / 1e4
             lamflamzp = 3943. * 3e-9 / 0.4206^2 * 0.4206	; Spitzer cookbook website (Fnu to Flam)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'VT': begin
             lameff(i)=5243.9 / 1e4
             weff(i)= 1001.6 / 1e4
             lamflamzp = 3761. * 3e-9 / 0.5244^2 * 0.5244	; Spitzer cookbook website (Fnu to Flam)
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'UKIS': begin
             lameff(i)=3581. / 1e4
             weff(i)= 638. / 1e4
             lamflamzp = 1.51e-5			; adopting same as U Johnson
             flux(i)=lamflamzp*10^(-0.4*mag(i))
           end
     'uSDSS': begin
             lameff(i)=poly(theta,[3472.,145.,77.,-55.]) / 1.e4
             weff(i)=poly(theta,[556.,263.,-562.,193.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.3530^2 * 0.3530
             flux(i)=lamflamzp*10^(-0.4*(mag(i)-0.04))  ; assuming Pogson mags instead of asinh
           end
     'gSDSS': begin
             lameff(i)=poly(theta,[4647.,312.,241.,-173.]) / 1.e4
             weff(i)=poly(theta,[1156.,909.,-1424.,387.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.4788^2 * 0.4788
             flux(i)=lamflamzp*10^(-0.4*mag(i))  ; assuming Pogson mags instead of asinh
           end
     'rSDSS': begin
             lameff(i)=poly(theta,[6145.,139.,156.,-80.]) / 1.e4
             weff(i)=poly(theta,[1255.,289.,-183.,-109.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.6242^2 * 0.6242
             flux(i)=lamflamzp*10^(-0.4*mag(i))  ; assuming Pogson mags instead of asinh
           end
     'iSDSS': begin
             lameff(i)=poly(theta,[7562.,101.,123.,-52.]) / 1.e4
             weff(i)=poly(theta,[1310.,144.,-9.,-104.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.7704^2 * 0.7704
             flux(i)=lamflamzp*10^(-0.4*mag(i))  ; assuming Pogson mags instead of asinh
           end
     'zSDSS': begin
             lameff(i)=poly(theta,[8997.,88.,105.,-36.]) / 1.e4
             weff(i)=poly(theta,[1357.,91.,24.,-76.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.9038^2 * 0.9038
             flux(i)=lamflamzp*10^(-0.4*(mag(i)+0.02))  ; assuming Pogson mags instead of asinh
           end
     'zPS' : begin
             lameff(i)=8660. / 1.e4
             weff(i)=720. / 1.e4
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             flux(i)=lamflamzp*10^(-0.4*(mag(i)))  ; assuming Pogson mags instead of asinh
           end
     'gKIS': begin
             lameff(i)=poly(theta,[4647.,312.,241.,-173.]) / 1.e4
             weff(i)=poly(theta,[1156.,909.,-1424.,387.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.4788^2 * 0.4788
             flux(i)=lamflamzp*10^(-0.4*(mag(i)-0.08))  ; assuming Pogson mags instead of asinh
           end
     'rKIS': begin
             lameff(i)=poly(theta,[6145.,139.,156.,-80.]) / 1.e4
             weff(i)=poly(theta,[1255.,289.,-183.,-109.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.6242^2 * 0.6242
             flux(i)=lamflamzp*10^(-0.4*(mag(i)+0.16))  ; assuming Pogson mags instead of asinh
           end
     'iKIS': begin
             lameff(i)=poly(theta,[7562.,101.,123.,-52.]) / 1.e4
             weff(i)=poly(theta,[1310.,144.,-9.,-104.]) / 1.e4
             lamflamzp = 3631. * 3e-9 / 0.7704^2 * 0.7704
             flux(i)=lamflamzp*10^(-0.4*(mag(i)+0.37))  ; assuming Pogson mags instead of asinh
           end
     'galNUV': begin
             lameff(i)=2267. / 1.e4
             weff(i)=732. / 1.e4
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
; http://galexgi.gsfc.nasa.gov/docs/galex/Documents/ERO_data_description_2.htm#_Toc58822546
           end
     'galFUV': begin
             lameff(i)=1596. / 1.e4    ; adjusted from 1516 A based on test with HD209458 (see papers/kelt3 directory)
             weff(i)=268. / 1.e4
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
           end
     'WISE1': begin
             lameff(i)=33526. / 1.e4
             weff(i)=6626. / 1.e4
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             if keyword_set(wiseab) then begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
             endif else begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)+2.683d0))  
             endelse
           end
     'WISE2': begin
             lameff(i)=46028. / 1.e4
             weff(i)=10423. / 1.e4
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             if keyword_set(wiseab) then begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
             endif else begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)+3.319d0))  
             endelse
           end
     'WISE3': begin
             lameff(i)=11.5608
             weff(i)=5.5069
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             if keyword_set(wiseab) then begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
             endif else begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)+5.242d0))  
             endelse
           end
     'WISE4': begin
             lameff(i)=22.0883
             weff(i)=4.1013
             lamflamzp = 3631. * 3e-9 / lameff(i)^2 * lameff(i)
             if keyword_set(wiseab) then begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)))  
             endif else begin
                flux(i)=lamflamzp*10^(-0.4*(mag(i)+6.604d0))  
             endelse
          end

;; JDE: 2018-07-13 -- I don't trust these...
     ;; Table 1 of Jordi et al, 2010
     ;; with DR2 Vega to AB mag offsets from
     ;; https://www.cosmos.esa.int/web/gaia/iow_20180316
     ;; commented offsets are more accurate, but not used for DR2 mags
;     'Gaia': begin ; Gaia
;             lameff(i)=0.673d0
;             weff(i)=0.440d0
;             lamflamzp = 3631d0 * 3d-9 / lameff(i)^2 * lameff(i)
;             flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.1050312311d0))  
;             ;flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.1001113078d0))  
;          end
;     'GBP': begin ; Gaia Blue
;             lameff(i)=0.532d0
;             weff(i)=0.253d0
;             lamflamzp = 3631d0 * 3d-9 / lameff(i)^2 * lameff(i)
;             flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.0291714680d0))  
;             ;flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.0373453185d0))  
;          end
;     'GRP': begin ; Gaia Red
;             lameff(i)=0.797d0
;             weff(i)=0.296d0
;             lamflamzp = 3631d0 * 3d-9 / lameff(i)^2 * lameff(i)
;             flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.3542076819d0))  
;             ;flux(i)=lamflamzp*10^(-0.4d0*(mag(i)+0.3534919681d0))  
;          end
     else: print, 'Band not recognized: ' + band[i]

  endcase
endfor

fluxerr = flux * alog(10)/2.5*merr

good = where(fluxerr ne 0)
if good[0] eq -1 then message, 'No good bands in SED file'

band = band[good]
lameff = lameff[good]
weff = weff[good]
flux = flux[good]
fluxerr = fluxerr[good]

;forprint,band,lameff,weff,flux,fluxerr

end

