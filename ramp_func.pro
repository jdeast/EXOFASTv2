pro ramp_func, x, a, f, pder

  ;; this is just a garbage function that should strongly penalize
  ;; negative exponential values (LC ramps don't behave that
  ;; way) in a way that should guide the fit toward negative values
  if a[1] ge 0d0 then begin
     f = 1d0+x*(10d0*(a[1]+1d0))
     if n_params() ge 4 then pder = [[x*0+1d0],[10d0*x]]
  endif else begin
     bx = exp(a[1]*x)
     f = 1d0+a[0]*bx
     if n_params() ge 4 then pder = [[bx],[a[0]*x*bx]]
  endelse  

end
