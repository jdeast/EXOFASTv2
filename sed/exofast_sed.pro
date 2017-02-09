function exofast_sed,fluxfile,teff,rstar,Av,d,logg=logg,met=met,f0=f0, fp0=fp0, ep0=ep0, verbose=verbose

common sed_block, klam, kkap, kapv, fp, ep, wp, widthhm, n, m

; fluxfile is the input set of fluxes to be fit
; a is the vector of input parameters
; f is the output model fluxes integrated over the bandpasses
; chi2 is output chi2

if not keyword_set(logg) then logg=4.5
if not keyword_set(met) then met=0.0
if not keyword_set(xmin) then xmin=0.1
if not keyword_set(xmax) then xmax=35.

c=2.9979e14
pc=3.0857e18
rsun=6.96e10

;; don't do this every time
if n_elements(fp) eq 0 then begin
   mag2fluxconv,fluxfile,wp,widthhm,fp,ep,teff=teff
   n=n_elements(wp)
   m=where(ep gt 0)
   readcol,getenv('EXOFAST_PATH') + '/sed/extinction_law.ascii',klam,kkap,/silent
   kapv = interpol(kkap,klam,0.55)
endif

;; round to 0.5 dex in each (!??!)
logg1 = double(round(logg*2d0))/2d0
met1 = double(round(met*2d0))/2d0

exofast_interp_model,teff,logg1,met1,w1,lamflam1temp,/next
if ~finite(lamflam1temp[0]) then return, !values.d_infinity
lamflam1=lamflam1temp*rstar*rstar*rsun*rsun/d/d/pc/pc
kapp1 = interpol(kkap,klam,w1)
taul1 = kapp1/kapv/1.086*Av
extinct1 = exp(-taul1)

flux = lamflam1*extinct1

f=fltarr(n)
for i=0,n-1 do begin
   ;; find the indices of the FWHM
   dummy = min((abs(w1-(wp[i]-widthhm[i]/2))),l)
   dummy = min((abs(w1-(wp[i]+widthhm[i]/2))),u)   
   f[i]=mean(flux[l:u])
endfor

if keyword_set(verbose) then begin
   if keyword_set(ps) then begin
      loadct, 39, /silent
      colors=[0,254,68,128]
   endif else begin
      device,window_state=win_state
      if win_state[5] eq 1 then wset, 5 $
      else window, 5
      colors = ['ffffff'x,'0000ff'x,'00ff00'x,'ff0000'x]
   endelse
   plotsym, 0, /fill, color=colors[1]
   plot, w1, flux,/xlog,/ylog,xtitle=textoidl('\lambda (\mum)'), ytitle=textoidl('log \lambda F_\lambda (erg s^{-1} cm^{-2})')
   oploterr, wp, f, ep, 8

endif

f0 = f
fp0 = fp
ep0 = ep
;chi2 = total( (alog10(f[m]) - alog10(fp[m]))^2 / (ep[m]/fp[m]*alog(10)/2.5)^2)
chi2 = total(((f-fp)/ep)^2)

return, chi2

end

