;+
; NAME:
;   CALCB
;
; PURPOSE: 
;   Calculates the impact parameter of a given time. The function to minimize
;
;   Implementation translated from here:
;   https://en.wikipedia.org/wiki/Golden-section_search
;
; CALLING SEQUENCE:
;    min = goldenratio('func', min, max)
;
; INPUTS:
;
;    FUNC - A string specifying the name of the function to be
;           minimized
;    A    - The minimum bound on the function to be minimized
;    B    - The maximum bound on the function to be minimized
;
; RESULT:
;
;    The minimum value
;
; MODIFICATION HISTORY
; 
;  2018/11/13 -- Written by Jason Eastman, CfA
;
;-
function calcb, time, logname=logname

  common tranblock, inc, ar, tp, period, e, omega, q

  return, exofast_getb2(time, inc=inc, a=ar, tperiastron=tp, period=period, $
                        e=e, omega=omega, q=q)

end

;+
; NAME:
;   EXOFAST_GETTT
; PURPOSE:
;   Calculates the time of minimum projected separation between the
;   star and planet (i.e., the time of transit).
;
; INPUTS:
;   MCMCSS - The MCMC stellar structure saved in the idl (*.mcmc.idl)
;            file output by exofastv2.pro  
;
; RESULT - An NPLANETS pointer array, pointing to an NEPOCHS x NSTEPS
;          array of Transit times, including stellar reflex motion and
;          fitted TTVs.
; 
;-
function exofast_gettt, mcmcss, epochs, secondary=secondary, tol=tol, filename=filename, filebase=filebase,debug=debug,nthin=nthin

common tranblock, inc, ar, tp, period, e, omega, q

if n_elements(filename) ne 0 then restore, filename
if n_elements(nthin) eq 0 then nthin = 1

if n_elements(tol) eq 0 then tol = 0.1d0/86400d0 ;; accurate to 0.1 seconds

chi2 = reform((*(mcmcss.chi2)),mcmcss.nsteps/mcmcss.nchains,mcmcss.nchains)
if (size(mcmcss.goodchains))[0] ne 10 then begin
   burnndx = getburnndx(chi2,goodchains=goodchains)
endif else begin
   burnndx = mcmcss.burnndx
   goodchains = (*mcmcss.goodchains)
endelse

;; nthin has to be a factor of the number of steps
nstepsperchain = mcmcss.nsteps/mcmcss.nchains
bestdiff = !values.d_infinity
for i=1L, nstepsperchain do begin
   if double(long(nstepsperchain/i)) eq double(nstepsperchain/i) then begin
      if abs(i-nthin) lt bestdiff then begin
         bestnthin = i
         bestdiff = abs(i-nthin)
      endif
   endif
endfor
nthin = bestnthin

if keyword_set(debug) then begin
   nsteps=mcmcss.nchains*100d0
   burnndx = 0
endif else nsteps = mcmcss.nsteps

;; **********************************************************************
;PRINT, 'this goes much faster for debugging by only computing a fraction of'
;PRINT, 'the steps -- DO NOT LEAVE THIS IN!!!'
;nsteps = mcmcss.nchains*5L ;; must be an integer multiple of nchains
;burnndx = 0
;; **********************************************************************

ttn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
bn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
depthn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
pn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
epochs = lonarr(mcmcss.nplanets,mcmcss.ntran)

;; secondary
tsn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
bsn = dblarr(nsteps,mcmcss.nplanets,mcmcss.ntran) + !values.d_nan
epochs2 = lonarr(mcmcss.nplanets,mcmcss.ntran)

for i=0L, nsteps-1 do begin

   ;; loop over each transit
   ;; each transit could (in principle) have information about every planet
   for k=0L, mcmcss.ntran-1 do begin
      
      ;; find the date span of this transit file
      minbjd0 = min((*(mcmcss.transit[k].transitptrs)).bjd,max=maxbjd0)

      ;; loop over each planet
      for j=0L, mcmcss.nplanets-1 do begin

         minbjd = minbjd0
         maxbjd = maxbjd0

         minbjds = minbjd0
         maxbjds = maxbjd0

         tp = mcmcss.planet[j].tp.value[i*nthin] + mcmcss.transit[k].ttv.value[i*nthin]      
         period = mcmcss.planet[j].period.value[i*nthin]

;         ;; shift to optimal epoch (shouldn't matter)
;         nper = round((mcmcss.planet[j].t0.value[i*nthin] - mcmcss.planet[j].tc.value[i*nthin])/period)
;         tp -= nper*period

         e = mcmcss.planet[j].e.value[i*nthin]
         omega = mcmcss.planet[j].omega.value[i*nthin]
         t14 = mcmcss.planet[j].t14.value[i*nthin]
         t14s = mcmcss.planet[j].t14s.value[i*nthin]
         inc = mcmcss.planet[j].i.value[i*nthin] + mcmcss.transit[k].tiv.value[i*nthin]
         ar = mcmcss.planet[j].ar.value[i*nthin]
         q = mcmcss.planet[j].q.value[i*nthin]
         p = mcmcss.planet[j].p.value[i*nthin]+mcmcss.transit[k].tdeltav.value[i*nthin]

         ;; primary transit
         phase = exofast_getphase(e,omega,/pri)
         tc = tp + period*phase
         epoch = round(((maxbjd+minbjd)/2d0 - tc)/period)
         thistc = tc + epoch*period

         ;; secondary occultation
         phase2 = exofast_getphase(e,omega,/sec)
         ts = tp + period*phase2
         epoch2 = round(((maxbjd+minbjd)/2d0 - ts)/period)
         thists = ts + epoch2*period

         ;; take the epoch from the middle of the run 
         ;; ***not robust if there is confusion about the period/epoch***
         if i eq nsteps/2 then begin
            epochs[j,k] = epoch
            epochs2[j,k] = epoch2
         endif

         ;; calculate T_T if the transit data overlaps the transit
         if maxbjd gt (thistc-t14/2d0) and minbjd lt (thistc+t14/2d0) then begin
            ;; find the minimum projected separation during the transit
            ;; assumes time of conjunction is the time of transit to
            ;; within the duration of the transit 
            ;; *** Might not be true if large bodies are in the system ***
            minbjd = (tc + period*epoch - t14)
            maxbjd = (tc + period*epoch + t14)
            ttn[i,j,k] = goldenratio('calcb',minbjd,maxbjd)
            
            ;; get the minimum projected separation (impact parameter)
            bn[i,j,k] = calcb(ttn[i,j,k])
            
            ;; get the limb darkened depth at the time of minimum
            ;; projected separation (T_T)
            u1 = mcmcss.band[mcmcss.transit[k].bandndx].u1.value[i*nthin]
            u2 = mcmcss.band[mcmcss.transit[k].bandndx].u2.value[i*nthin]
            exofast_occultquad_cel, bn[i,j,k], u1, u2, p, mu1
            depthn[i,j,k] = 1d0-mu1
            pn[i,j,k] = p

         endif

         ;; calculate T_{S,T} if the occultation data overlaps the transit
         if maxbjds gt (thists-t14s/2d0) and minbjds lt (thists-t14s/2d0) then begin
        
            ;; find the minimum projected separation during the transit
            ;; assumes time of conjunction is the time of transit to
            ;; within the duration of the transit 
            ;; *** Might not be true if large bodies are in the system ***
            minbjds = (ts + period*epoch2 - t14s)
            maxbjds = (ts + period*epoch2 + t14s)
            tsn[i,j,k] = goldenratio('calcb',minbjds,maxbjds)
            
            ;; get the minimum projected separation (impact parameter)
            bsn[i,j,k] = calcb(tsn[i,j,k])
            
         endif

      endfor
   endfor

   if ((i+1d0) mod long(nsteps/10d0)) eq 0 then printandlog, strtrim(round((i+1d0)/nsteps*100,/L64),2) + '% done computing transit times', logname

endfor

;; nsteps x nchains x nplanets x ntransits array
ttn = reform(ttn,nsteps/mcmcss.nchains,mcmcss.nchains,mcmcss.nplanets,mcmcss.ntran) 
bn = reform(bn,nsteps/mcmcss.nchains,mcmcss.nchains,mcmcss.nplanets,mcmcss.ntran) 
depthn = reform(depthn,nsteps/mcmcss.nchains,mcmcss.nchains,mcmcss.nplanets,mcmcss.ntran) 
pn = reform(pn,nsteps/mcmcss.nchains,mcmcss.nchains,mcmcss.nplanets,mcmcss.ntran) 
medianttn = [-1]
medianbn = [-1]
mediandepthn = [-1]
medianpn = [-1]
planetname = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

if n_elements(filebase) eq 0 then filebase = 'planet.'
openw, texlun, filebase + 'transits.tex', /get_lun
openw, csvlun, filebase + 'transits.csv', /get_lun
printf, csvlun, '# Transit, Planet, Epoch, T_T, sigma_tt_high, sigma_tt_low, b, sigma_b_high, sigma_b_low, Depth, sigma_depth_high, sigma_depth_low, p, sigma_p_high, sigma_p_low'

printf, texlun, '\documentclass{aastex62}'
printf, texlun, '\begin{document}'
printf, texlun, '\startlongtable'
printf, texlun, '\begin{deluxetable*}{lcccccc}'
printf, texlun, '\tablecaption{Median values and 68\% confidence interval for transit times, impact parameters, and depths}'
printf, texlun, '\tablehead{\colhead{Transit} & \colhead{Planet} & \colhead{Epoch} & \colhead{$T_T$} & \colhead{$b$} & \colhead{Depth} & \colhead{$R_P/R_*$}}'
printf, texlun, '\startdata'

for j=0L, mcmcss.nplanets-1 do begin
   for k=0L, mcmcss.ntran-1 do begin
      good = where(finite(ttn[burnndx:*,goodchains,j,k]),ngood)
      if ngood gt 0 then begin
         summarizepar, ttn[burnndx:*,goodchains,j,k], medianpars=medianttn, value=thistt, errlo=ttlo, errhi=tthi, logname=logname, /noplot
         summarizepar, bn[burnndx:*,goodchains,j,k], medianpars=medianbn, value=thisb, errlo=blo, errhi=bhi, logname=logname, /noplot
         summarizepar, depthn[burnndx:*,goodchains,j,k], medianpars=mediandepthn, value=thisdepth, errlo=depthlo, errhi=depthhi, logname=logname, /noplot
         summarizepar, pn[burnndx:*,goodchains,j,k], medianpars=medianpn, value=thisp, errlo=plo, errhi=phi, logname=logname, /noplot
         printf, csvlun, mcmcss.transit[k].label, planetname[j], strtrim(epochs[j,k],2), strtrim(thistt,2), strtrim(tthi,2), strtrim(ttlo,2), strtrim(thisb,2), strtrim(bhi,2), strtrim(blo,2), strtrim(thisdepth,2), strtrim(depthhi,2), strtrim(depthlo,2), strtrim(thisp,2), strtrim(phi,2), strtrim(plo,2), format='(a,",",a1,",",a,12(",",a))'

         ;; if the error bars are the same, output +/- instead of ^+_-
         if tthi ne ttlo then tt = string(strtrim(thistt,2), strtrim(tthi,2), strtrim(ttlo,2), format='("$",a,"^{+",a,"}_{-",a,"}$")') $
         else tt = string(strtrim(thistt,2), strtrim(tthi,2), format='("$",a," \pm ",a,"$")')
         if bhi ne blo then b = string(strtrim(thisb,2), strtrim(bhi,2), strtrim(blo,2), format='("$",a,"^{+",a,"}_{-",a,"}$")') $
         else b = string(strtrim(thisb,2), strtrim(bhi,2), format='("$",a," \pm ",a,"$")')
         if depthhi ne depthlo then depth = string(strtrim(thisdepth,2), strtrim(depthhi,2), strtrim(depthlo,2), format='("$",a,"^{+",a,"}_{-",a,"}$")') $
         else depth = string(strtrim(thisdepth,2), strtrim(depthhi,2), format='("$",a," \pm ",a,"$")')
         if phi ne plo then p = string(strtrim(thisp,2), strtrim(phi,2), strtrim(plo,2), format='("$",a,"^{+",a,"}_{-",a,"}$")') $
         else p = string(strtrim(thisp,2), strtrim(phi,2), format='("$",a," \pm ",a,"$")')

         printf, texlun, mcmcss.transit[k].label, planetname[j], strtrim(epochs[j,k],2), tt, b, depth, p, format='(a," & ",a1," & ",a,4(" & ",a),"\\")'
      endif
   endfor
endfor

printf, texlun, '\enddata'
printf, texlun, '\label{tab:transitpars}'
printf, texlun, '\end{deluxetable*}'
printf, texlun, '\end{document}'
free_lun, texlun
free_lun, csvlun

return, ttn

end

