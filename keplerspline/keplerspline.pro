function keplerspline2, t, f, ndays=ndays, maxiter = maxiter, rms = rms, yplot = yplot, goodout = goodout

  t2 = (t - min(t)) / (max(t) - min(t))
  bksp = ndays / (max(t) - min(t))
  
  lastgood = findgen(n_elements(t2))
  
  for i = 0, maxiter do begin
  
  
    thist = t2(lastgood)
    thisf = f(lastgood)
    
    sset = bspline_iterfit(thist,thisf,nord=4,maxiter=0,bkspace=bksp)
    
    bg = bspline_valu(t2,sset)
    
;    if total(1 - finite(bg)) gt 0 then begin
;      stop
;      return, dblarr(n_elements(bg)) + 1d
;    endif
    
    junk = robust_mean(f - bg, 3, goodind = good)
    rms = stddev(bg(good) - f(good))
    if keyword_set(yplot) then begin
    
      if i le yplot then begin
        small = 0.7
        large = 1
        title = 'Iteration '+ trim(i + 1, 1)
        if i eq 0 then title = 'EPIC 60019950: Iteration '+ trim(i + 1, 1)
        if i eq 2 then xtitle = 'BJD - 2454833'
        
        plot, t, f, ps=8, /ynozero,/nodata, title = title,ytitle = 'Relative Flux', xtitle = xtitle
        oplot, t, f, ps=8, symsize = large
        oplot, t, f, ps=8, symsize = small, color = cgcolor('red')
        oplot, t(lastgood), f(lastgood), ps=8;, color = cgcolor('black'), symsize = large
        oplot, t, bg, thick = 15;, color = cgcolor('black')
        oplot, t, bg, thick = 12, color = cgcolor('orange')
      endif
    endif
    if n_elements(lastgood) eq n_elements(good) then begin
      if total(lastgood ne good) eq 0 then begin
        if keyword_set(yplot) then print, 'converged in ', i + 1 , ' iterations'
        break
      endif
    end
    
    
    lastgood = good
    
  end
  goodout = good
  return, bg
  
end
function keplerspline, t, f, ndays=ndays, maxiter = maxiter, rms = rms, yplot = yplot, breakp = breakp, goodind = goodind


  if 1 - keyword_set(ndays) then ndays = 1.5
  if n_elements(maxiter) eq 0 then maxiter = 5

  gaps = where(diff(t) gt ndays)
  if keyword_set(breakp) then begin
    if gaps[0] eq -1 then gaps = [breakp]
    if gaps[0] gt -1 then begin
      gaps = [gaps, breakp]
      gaps = gaps[UNIQ(gaps, SORT(gaps))]
    endif
  endif

  if gaps[0] eq -1 then begin
    s = keplerspline2(t, f, ndays = ndays, maxiter = maxiter, rms=rms, yplot=yplot, good = good) ; make sure we separate the continuum normalizations between gaps in data
    goodind = good
  endif
  if gaps[0] gt -1 then begin
    s = dblarr(n_elements(t))
    goodind = [-1]
    gaps = [0,gaps,n_elements(t)]
    for i =1, n_elements(gaps)-1 do begin
      ;print, max(diff(k.t[gaps[i-1]:gaps[i]-1]))
      s[gaps[i-1]:gaps[i]-1] = keplerspline2(t[gaps[i-1]:gaps[i]-1], f[gaps[i-1]:gaps[i]-1], ndays = ndays, maxiter = maxiter,rms=rms,yplot=yplot, good = good)
      goodind = [goodind, good + gaps[i-1]]
    endfor
    if n_elements(goodind) gt 1 then goodind = goodind[1:*]
  end
  
  
  return, s
end

