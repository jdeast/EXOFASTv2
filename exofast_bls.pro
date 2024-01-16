function exofast_bls, time, flux, minperiod=minperiod, maxperiod=maxperiod, $
                      nperiod=nperiod, nbins=nbins, maxfrac=maxfrac, plot=plot, $
                      powers=powers,t10=t10,t4=t4,t14=t14,depth=depth,ndx=ndx,periods=periods, tc=tc

nflux = n_elements(flux)

if n_elements(minperiod) eq 0 then minperiod = 0.5d
if n_elements(maxperiod) eq 0 then maxperiod = 1d1
if n_elements(nperiod) eq 0 then nperiod = 3d4
if n_elements(nbins) eq 0 then nbins = 300d
if n_elements(maxfrac) eq 0 then maxfrac = 0.25d0

;; look for periods evenly spaced in frequency
minfreq = 1d0/maxperiod
maxfreq = 1d0/minperiod
freq = minfreq + (maxfreq-minfreq)*dindgen(nperiod)/(nperiod-1d0)
periods = reverse(1d0/freq)

;; can't look for periods longer than the range
timespan = time[nflux-1]-time[0]
if maxperiod gt timespan then message, 'Period > Timespan'

maxdur = maxfrac*nbins + 1d0

;; the average time
tave = mean(time) 

zeroaveflux = flux-mean(flux) 
powers = dblarr(nperiod)
t10 = dblarr(nperiod)
t4 = dblarr(nperiod)
t14 = dblarr(nperiod)
depth = dblarr(nperiod)

for i=0L, nperiod-1 do begin

    ;; the phase of each time
    phase = (time mod periods[i])/periods[i]
    
    ;; get the total zero-averaged flux in each bin
    binflux = dblarr(nbins)
    hist = histogram(phase,nbins=nbins,min=0,max=(1d0-1d0/nbins),$
                     reverse_ind=R,locations=binphase)
    for j=0L, nbins-1 do $
      if R[j] ne R[j+1] then binflux[j] = total(zeroaveflux[R[R[j]:R[j+1]-1]])
    
    ;; extend the arrays to compare the boundaries
    hist = [hist,hist[0:maxdur]]
    binflux = [binflux,binflux[0:maxdur]]

    ;; for each phase (starting point)
    for j=0L, nbins-1 do begin

        ;; calculate the power at each trial duration
        s = total(binflux[j:j+maxdur-1],/cumulative)
        kk = total(hist[j:j+maxdur-1],/cumulative)
        pow = s^2/(nflux*(nflux-kk))
        
        ;; if it's the best bin for this period, save it
        maxpow = max(pow,ndx)
        if maxpow gt powers[i] then begin

            ;; calculate interesting bits about this power 
            ;; (transit duration, tc, depth)
            t0 = tave - (tave mod periods[i])

            t10[i] = t0 + binphase[j]*periods[i]
            t4[i] = t0 + binphase[(j+ndx) mod nbins]*periods[i]
            t14[i] = kk[ndx]/nflux*periods[i]
            depth[i] = -s[ndx]*nflux/(kk[ndx]*(nflux-kk[ndx]))
            powers[i] = maxpow

        endif

    endfor

endfor

;; the power is actually the sqrt
powers = sqrt(powers)

;; central transit times
tc = (t4+t10)/2d0
toohigh = where(t10 gt t4)
if toohigh[0] ne -1 then tc[toohigh] += periods[toohigh]/2d0

;; best period
maxpow = max(powers,ndx)
bestperiod = periods[ndx]

if keyword_set(plot) then begin

    !p.multi = [0,2]
    phase = (((time - tc[ndx]) mod periods[ndx])/periods[ndx] + 1.5d0) mod 1
    plot, phase, flux, psym=1, /ystyle
    oplot, [0.5,0.5],[-9d9,9d9], color='0000ff'x
    oplot, [-9d9,9d9],[1,1]
    oplot, [-9d9,9d9],[1-depth[ndx],1-depth[ndx]]

    plot, periods, powers;, /xlog, /xstyle
    !p.multi=0

endif

return, bestperiod

end



